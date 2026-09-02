#!/usr/bin/env python3
import os
import sys
import time
import pgenlib
import tempfile
import warnings
import numpy as np
from pathlib import Path
import concurrent.futures
import multiprocessing as mpi
from test_pgenlib import unphased_biallelic_case

try:
    num_cpus = len(os.sched_getaffinity(os.getpid()))
except AttributeError:
    # if on macos, fallback to number of CPUs given by os.cpu_count()
    num_cpus = mpi.cpu_count()
print(f"Using {num_cpus} CPUs", file=sys.stderr)

def generate_large_pgen(pgen_dir, case_idx=0, nsample_min=1, nsample_limit=20000, nvariant_min=1, nvariant_limit=80000):
    pgen_dir = Path(pgen_dir)
    if pgen_dir.exists():
        return next(pgen_dir.glob("*.pgen"))  # Return existing .pgen file if present
    pgen_dir.mkdir(parents=True, exist_ok=True)  # Ensure directory exists
    # Generate a single large pgen file
    start = time.time()
    unphased_biallelic_case(pgen_dir, case_idx, nsample_min, nsample_limit, nvariant_min, nvariant_limit)
    print(f"Time to generate PGEN: {time.time() - start:.2f}s", file=sys.stderr)
    # Return the expected .pgen file path
    return next(pgen_dir.glob("*.pgen"))

def timed_read(file_path, return_arr=True):
    with pgenlib.PgenReader(str(file_path).encode()) as r:
        num_vars = r.get_variant_ct()
        arr = np.empty([num_vars, r.get_raw_sample_ct()*2], dtype=np.int32)
        start = time.time()
        r.read_alleles_range(0, num_vars, arr)
        elapsed = time.time() - start
    if return_arr:
        return elapsed, arr
    else:
        return elapsed

def threaded_timed_read(file_path, n_threads=num_cpus, single_arr=None):
    file_path = str(file_path).encode()
    with pgenlib.PgenReader(file_path) as r:
        num_vars, n_samples = r.get_variant_ct(), r.get_raw_sample_ct()
        arr = np.empty((num_vars, n_samples * 2), dtype=np.int32)

    def do_reads(s, e):
        with pgenlib.PgenReader(file_path) as o:
            o.read_alleles_range(s, e, arr[s:e])

    start = time.time()
    chunk = (num_vars + n_threads - 1) // n_threads
    with concurrent.futures.ThreadPoolExecutor(n_threads) as ex:
        ex.map(lambda i: do_reads(i * chunk, min((i + 1) * chunk, num_vars)), range(n_threads))
    elapsed = time.time() - start
    if single_arr is not None:
        np.testing.assert_allclose(single_arr, arr)
    return elapsed

def process_init_mp(
    file_path_: str = None,
    num_vars_: int = None,
    n_samples_: int = None,
    shared_arr_: mpi.Array = None,
):
    """
    A helper method that globalizes certain variables so that they can be used in
    multiprocessing. This method should only be called in each parallel
    child/worker process and not the parent process to avoid polluting the global
    namespace.

    This is preferable to passing these values as arguments to the _do_reads method
    because, otherwise, python will try to pickle the arguments

    Parameters
    ----------
    shared_arr: mpi.Array
        The underlying bytes of the matrix, as a shared-memory Array
    """
    global shared_arr, num_vars, n_samples, file_path
    shared_arr, num_vars, n_samples, file_path = shared_arr_, num_vars_, n_samples_, file_path_

def process_do_reads(s, e):
    shd_arr = shared_arr.get_obj()
    with pgenlib.PgenReader(file_path) as o:
        arr = np.frombuffer(shd_arr, dtype=np.int32).reshape((num_vars, n_samples * 2))
        o.read_alleles_range(s, e, arr[s:e])

def process_timed_read(file_path, n_threads=num_cpus, single_arr=None):
    return np.inf # comment this out to enable multi-process testing
    file_path = str(file_path).encode()
    with pgenlib.PgenReader(file_path) as r:
        num_vars, n_samples = r.get_variant_ct(), r.get_raw_sample_ct()
        shared_arr = mpi.Array("i", int(np.prod((num_vars, n_samples * 2))))
        shd_arr = shared_arr.get_obj()
        arr = np.frombuffer(shd_arr, dtype=np.int32).reshape((num_vars, n_samples * 2))

    start = time.time()
    chunk = (num_vars + n_threads - 1) // n_threads
    pairs = [(i * chunk, min((i + 1) * chunk, num_vars)) for i in range(n_threads)]
    mp_chunksize = int(np.ceil(len(pairs) / n_threads))
    with mpi.Pool(
        processes=n_threads,
        initargs=(file_path, num_vars, n_samples, shared_arr),
        initializer=process_init_mp,
    ) as ex:
        ex.starmap(process_do_reads, pairs, chunksize=mp_chunksize)
    elapsed = time.time() - start
    if single_arr is not None:
        np.testing.assert_allclose(single_arr, arr)
    return elapsed

def bench_instance(pgen_path: Path):
    single, single_arr = timed_read(pgen_path)
    thread = threaded_timed_read(pgen_path, single_arr=single_arr)
    proc = process_timed_read(pgen_path, single_arr=single_arr)

    return single, thread, proc, single_arr.shape

def ols_fit_1d(x: np.ndarray, y: np.ndarray, reps: int = 1):
    """
    Fit y = intercept + slope * x by ordinary least squares using numpy.
    Returns (intercept, slope).
    """
    if y.ndim > 1:
        x = np.repeat(np.asarray(x, dtype=float).ravel(), reps)
    else:
        x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    if x.size != y.size:
        raise ValueError("x and y must have the same length")
    A = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    intercept, slope = coef[0], coef[1]

    return float(intercept), float(slope)

def main(tmp_path, also_plot: bool = True):
    nsample_fixed = 200000
    nvariant_limits = np.linspace(5, 118, 5, dtype=int)
    num_reps = 7
    single_times = []
    thread_times = []
    process_times = []
    shapes = []
    valid_result = True

    nsample_actual = None
    for nvariant_limit in nvariant_limits:
        pgen_dir = tmp_path / f"sample_{nvariant_limit}"
        # generate_large_pgen is assumed to return a path to the generated pgen file
        pgen_path = generate_large_pgen(pgen_dir, nvariant_limit=nvariant_limit, nsample_limit=nsample_fixed)

        single_times_rep = []
        thread_times_rep = []
        process_times_rep = []
        shapes_rep = set()

        for rep in range(num_reps):
            s_time, t_time, p_time, shape = bench_instance(pgen_path)
            single_times_rep.append(s_time)
            thread_times_rep.append(t_time)
            process_times_rep.append(p_time)
            shapes_rep.add(shape)

        s_mean = float(np.mean(single_times_rep))
        t_mean = float(np.mean(thread_times_rep))
        p_mean = float(np.mean(process_times_rep))
        s_med = float(np.median(single_times_rep))
        t_med = float(np.median(thread_times_rep))
        p_med = float(np.median(process_times_rep))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            s_std = float(np.std(single_times_rep))
            t_std = float(np.std(thread_times_rep))
            p_std = float(np.std(process_times_rep))

        assert len(shapes_rep) == 1, "Inconsistent shapes across repetitions"
        nvars, nsamps = shapes_rep.pop()
        if nsample_actual is None:
            nsample_actual = int(nsamps // 2)
        else:
            assert nsample_actual == int(nsamps // 2), "Inconsistent sample counts across different nvariant runs"

        single_times.append((s_mean, s_std, s_med))
        thread_times.append((t_mean, t_std, t_med))
        process_times.append((p_mean, p_std, p_med))
        shapes.append(nvars)

        # invalidate result if using mean or median would give different conclusions
        invalid = ((s_mean-t_mean)/abs(s_mean-t_mean)) != ((s_med-t_med)/abs(s_med-t_med))
        if s_std > 0.001 and invalid:
            valid_result = False
            print(
                f"nvariants={nvars}: ",
                f"single={s_med:.3f} {s_mean:.3f} ± {s_std:.3f} ({sorted(single_times_rep)})",
                f"thread={t_med:.3f} {t_mean:.3f} ± {t_std:.3f} ({sorted(thread_times_rep)})",
                (f"process={p_med:.3f} {p_mean:.3f} ± {p_std:.3f} ({sorted(process_times_rep)})" if p_mean < np.inf else ""),
                file=sys.stderr,
            )
            print("Mean and median disagree. Invalidating result.", file=sys.stderr)
        else:
            print(
                f"nvariants={nvars}: ",
                f"single={s_med:.3f} {s_mean:.3f} ± {s_std:.3f}",
                f"thread={t_med:.3f} {t_mean:.3f} ± {t_std:.3f}",
                (f"process={p_med:.3f} {p_mean:.3f} ± {p_std:.3f}" if p_mean < np.inf else ""),
                file=sys.stderr,
            )

    # Convert to numpy arrays for numeric ops
    X = np.array(shapes)
    single_times = np.array(single_times, dtype=float)
    thread_times = np.array(thread_times, dtype=float)
    process_times = np.array(process_times, dtype=float)

    # Fit OLS models using numpy-only ols_fit_1d
    models = {
        "Single-threaded": (single_times, ols_fit_1d(X, single_times[:,2], reps=num_reps)),
        "Multi-process": (process_times, ols_fit_1d(X, process_times[:,2], reps=num_reps)),
        "Multi-threaded": (thread_times, ols_fit_1d(X, thread_times[:,2], reps=num_reps)),
    }

    if also_plot:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(8, 6))
        for label, (times, (intercept, slope)) in models.items():
            # skip degenerate slope near zero
            if np.isnan(slope):
                continue
            # Plot with points at median (times[:,2])
            line, = plt.plot(X, times[:,2], 'o', label=f"{label} (slope={slope:.5e})")
            # Plot error bars centered at mean (times[:,0])
            plt.errorbar(X, times[:,0], yerr=times[:,1], fmt='none', color=line.get_color(), alpha=0.5)
            y_pred = intercept + slope * X
            plt.plot(X, y_pred, '--', color=line.get_color())

        plt.xlabel(f"Number of variants\nNumber of samples fixed at {nsample_actual}")
        plt.ylabel('Read time (s)')
        plt.title(f"PgenReader Benchmark with {num_cpus} CPUs")
        plt.legend(loc='upper left')
        plt.tight_layout()
        output_file_path = str(tmp_path / 'timings_plot.png')
        plt.savefig(output_file_path)
        print(f"Plot saved as {output_file_path}", file=sys.stderr)

    # check that multi-threaded is faster than both single-threaded and multi-processing
    # by comparing the slopes of the lines
    thread_slope = float(models["Multi-threaded"][1][1])
    print(f"Multi-threaded slope: {thread_slope}", file=sys.stderr)
    for model in ("Single-threaded", "Multi-process"):
        slope = float(models[model][1][1])
        if not np.isnan(slope):
            print(f"{model} slope: {slope}", file=sys.stderr)
            # skip the test if we don't have at least 2 CPUs or
            # if the err was too high for one of the points
            if num_cpus > 1 and valid_result:
                assert (slope > thread_slope)

def test_multithread(tmp_path):
    main(tmp_path, also_plot=False)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(Path(sys.argv[1]))
    else:
        with tempfile.TemporaryDirectory() as temp_dir:
            main(Path(temp_dir))

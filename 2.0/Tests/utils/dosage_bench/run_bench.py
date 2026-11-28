#!/usr/bin/env python3
"""Dosage benchmark: times typical plink2 command lines on dosage data in
three configurations:
  master-bi:  a master plink2 binary on an all-biallelic fileset
  dev-bi:     a dev plink2 binary on the same fileset
  dev-multi:  the dev binary on the partly-multiallelic version of it
The filesets come from gen_dosage_data.cc (built here), which derives the
all-biallelic fileset from the partly-multiallelic one by collapsing every
non-REF allele into ALT1; see its header comment for the genotype model.

Example:
  python3 run_bench.py --master /path/to/master/plink2 \\
      --dev /path/to/dev/plink2 --work-dir /tmp/dosage_bench

Each command runs once untimed (warm file cache) and then --reps times; the
median wall time and the max peak RSS are reported, along with the
dev-bi/master-bi and dev-multi/dev-bi ratios.  Results also go to
<work dir>/results.tsv.  Generated data is reused when the sizes and
generator options match (see <work dir>/data/params.txt).
"""

import argparse
import os
import shlex
import statistics
import subprocess
import sys
import time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INCLUDE_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "..", "..", "include"))

# name -> extra arguments; {data} is the generated-data prefix (for the
# phenotype/covariate/score/keep files), {out} the output prefix.
COMMANDS = [
    ("freq", "--freq"),
    ("make-pgen-keep", "--keep {data}.keep --make-pgen"),
    ("export-vcf-ds", "--export vcf bgz vcf-dosage=DS"),
    ("glm-linear", "--pheno {data}.pheno --pheno-name QT --covar {data}.covar --glm hide-covar"),
    ("glm-logistic", "--pheno {data}.pheno --pheno-name CC --covar {data}.covar --glm hide-covar firth-fallback"),
    ("score", "--score {data}.score 1 2 3 header cols=+scoresums"),
    ("pca-approx", "--pca approx 10"),
]

CONFIGS = [("master-bi", "master", "bi"), ("dev-bi", "dev", "bi"), ("dev-multi", "dev", "multi")]


def build_generator(work_dir):
    exe = os.path.join(work_dir, "gen_dosage_data")
    src = os.path.join(SCRIPT_DIR, "gen_dosage_data.cc")
    srcs = [src] + [os.path.join(INCLUDE_DIR, f) for f in
                    ("pgenlib_write.cc", "pgenlib_misc.cc", "plink2_base.cc", "plink2_bits.cc", "plink2_simd.cc")]
    if os.path.exists(exe) and all(os.path.getmtime(exe) > os.path.getmtime(s) for s in srcs):
        return exe
    cxx = os.environ.get("CXX", "c++")
    cmd = [cxx, "-std=c++17", "-O2", "-I" + INCLUDE_DIR, "-o", exe] + srcs + ["-lz"]
    print("building generator:", " ".join(shlex.quote(c) for c in cmd), flush=True)
    subprocess.run(cmd, check=True)
    return exe


def generate(args, work_dir):
    data_dir = os.path.join(work_dir, "data")
    os.makedirs(data_dir, exist_ok=True)
    prefix = os.path.join(data_dir, "bench")
    params = "%d %d %s" % (args.samples, args.variants, " ".join(args.gen_opt))
    params_fname = os.path.join(data_dir, "params.txt")
    if os.path.exists(params_fname) and os.path.exists(prefix + "_bi.pgen"):
        with open(params_fname) as f:
            if f.read().strip() == params:
                print("reusing data in", data_dir, flush=True)
                return prefix
    exe = build_generator(work_dir)
    cmd = [exe, prefix, str(args.samples), str(args.variants)] + args.gen_opt
    print("generating:", " ".join(cmd), flush=True)
    start = time.time()
    subprocess.run(cmd, check=True)
    print("  %.1f s" % (time.time() - start), flush=True)
    with open(params_fname, "w") as f:
        f.write(params + "\n")
    return prefix


def run_once(cmd, log_fname):
    """Returns (wall seconds, peak RSS in MiB)."""
    with open(log_fname, "w") as log:
        start = time.perf_counter()
        proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
        _, status, rusage = os.wait4(proc.pid, 0)
        wall = time.perf_counter() - start
    proc.returncode = os.waitstatus_to_exitcode(status) if hasattr(os, "waitstatus_to_exitcode") else (status >> 8)
    if proc.returncode != 0:
        sys.exit("command failed (exit code %d; see %s): %s" % (proc.returncode, log_fname, " ".join(cmd)))
    # ru_maxrss is in bytes on macOS, KiB on Linux
    rss_mib = rusage.ru_maxrss / (1048576.0 if sys.platform == "darwin" else 1024.0)
    return wall, rss_mib


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--master", required=True, help="master plink2 binary")
    parser.add_argument("--dev", required=True, help="dev plink2 binary")
    parser.add_argument("--work-dir", required=True, help="directory for data, outputs, and logs")
    parser.add_argument("--samples", type=int, default=10000)
    parser.add_argument("--variants", type=int, default=20000)
    parser.add_argument("--gen-opt", action="append", default=[],
                        help="generator option, e.g. multi=0.3 (repeatable)")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--memory", type=int, default=8000, help="plink2 --memory value (MiB)")
    parser.add_argument("--reps", type=int, default=3)
    parser.add_argument("--commands", help="comma-separated subset of: " + ",".join(c[0] for c in COMMANDS))
    args = parser.parse_args()

    work_dir = os.path.abspath(args.work_dir)
    os.makedirs(work_dir, exist_ok=True)
    data = generate(args, work_dir)
    bins = {"master": os.path.abspath(args.master), "dev": os.path.abspath(args.dev)}
    commands = COMMANDS
    if args.commands:
        wanted = args.commands.split(",")
        unknown = set(wanted) - set(c[0] for c in COMMANDS)
        if unknown:
            sys.exit("unknown command(s): " + ",".join(sorted(unknown)))
        commands = [c for c in COMMANDS if c[0] in wanted]
    out_dir = os.path.join(work_dir, "out")
    os.makedirs(out_dir, exist_ok=True)

    results = []
    for name, extra in commands:
        row = {}
        for config, bin_key, dataset in CONFIGS:
            out = os.path.join(out_dir, "%s_%s" % (name, config))
            cmd = [bins[bin_key], "--pfile", "%s_%s" % (data, dataset), "--threads", str(args.threads),
                   "--memory", str(args.memory)] + shlex.split(extra.format(data=data)) + ["--out", out]
            log_fname = out + ".stdout"
            run_once(cmd, log_fname)
            walls = []
            peak_rss = 0.0
            for _ in range(args.reps):
                wall, rss = run_once(cmd, log_fname)
                walls.append(wall)
                peak_rss = max(peak_rss, rss)
            row[config] = (statistics.median(walls), peak_rss)
            print("%-15s %-10s %8.2f s  %7.0f MiB" % (name, config, row[config][0], peak_rss), flush=True)
        results.append((name, row))

    header = ["command"] + ["%s_s" % c[0] for c in CONFIGS] + ["dev/master", "multi/bi"] + ["%s_MiB" % c[0] for c in CONFIGS]
    lines = ["\t".join(header)]
    for name, row in results:
        t = [row[c[0]][0] for c in CONFIGS]
        m = [row[c[0]][1] for c in CONFIGS]
        lines.append("\t".join([name] + ["%.2f" % x for x in t] + ["%.2f" % (t[1] / t[0]), "%.2f" % (t[2] / t[1])] +
                               ["%.0f" % x for x in m]))
    with open(os.path.join(work_dir, "results.tsv"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print()
    widths = [max(len(line.split("\t")[i]) for line in lines) for i in range(len(header))]
    for line in lines:
        print("  ".join(field.rjust(widths[i]) if i else field.ljust(widths[i]) for i, field in enumerate(line.split("\t"))))


if __name__ == "__main__":
    main()

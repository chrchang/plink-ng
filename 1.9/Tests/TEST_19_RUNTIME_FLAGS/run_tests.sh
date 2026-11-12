#!/bin/bash

# PLINK 1.9 runtime and environment flags: --version, --silent, --debug,
# --script, --rerun and --perm-batch-size.
#
# Usage: ./run_tests.sh {plink 1.9 dir}
# plink (1.9) is taken from {plink 1.9 dir}, or from the system PATH if that
# is empty or omitted.  One --rerun check reads a PLINK 1.07 log; the 1.07
# binary is taken from $PLINK107, defaulting to ../../../1.07/plink, and that
# check is skipped if it isn't there.

set -uo pipefail

# plink 1.9 from {plink 1.9 dir}, or from the PATH if that is empty or omitted.
plink19="${1:+$1/}plink"
# (made absolute, since one check runs in a subdirectory)
if [[ "$plink19" == */* && "${plink19:0:1}" != "/" ]]; then
    plink19="$PWD/$plink19"
fi
plink() {
    command "$plink19" "$@"
}
plink107="${PLINK107:-../../../1.07/plink}"

failed=0
check() {
    if [[ "$2" -eq 0 ]]; then
        echo "ok: $1"
    else
        echo "FAILED: $1"
        failed=1
    fi
}

# A log without the lines that legitimately differ between two runs of the
# same command: the flag list, host, directory and times.  The output prefix
# ($2) is replaced with OUT.
log_body() {
    grep -v '^  --\|^Hostname:\|^Working directory:\|^Start time:\|^End time:' "$1" | sed "s/$2/OUT/g"
}

# The flag list ("Options in effect:" block) of a log.
log_flags() {
    awk '/^Options in effect:/ { on = 1; next } on && !/^  / { exit } on' "$1"
}

# The flag list of log $1, with output prefix $2 replaced by OUT and flag $3
# (if given) added, sorted.
flags_plus() {
    { log_flags "$1"; [[ -n "${3:-}" ]] && echo "  $3"; } | sed "s/$2/OUT/" | LC_ALL=C sort
}

# 240 unrelated samples over 400 variants on two chromosomes, with case/
# control and quantitative phenotypes; 60 families of four with quantitative
# phenotypes; and 8 sets of 50 consecutive variants.
awk -v n=240 -v m=400 -v seed=7 -f make_data.awk > tmp_cc.ped
awk -v n=240 -v m=400 -v seed=7 -v qt=1 -f make_data.awk > tmp_qt.ped
awk -v n=240 -v m=400 -v seed=9 -v qt=1 -v fam=1 -f make_data.awk > tmp_fq.ped
awk 'BEGIN { for (j = 1; j <= 400; ++j) print 1 + int((j - 1) / 200), "v" j, 0, 1000 * j }' > tmp_cc.map
cp tmp_cc.map tmp_qt.map
cp tmp_cc.map tmp_fq.map
for f in cc qt fq; do
    plink --file tmp_$f --make-bed --out tmp_$f > /dev/null 2>&1
done
awk '{ print $1, $4, $4, "S" int((NR - 1) / 50) }' tmp_cc.bim > tmp_sets.txt

# --- --version ---

# -v and -V are aliases.  Other flags are ignored: nothing is computed and no
# log is written.
ver_ok=0
plink --version > tmp_ver.txt 2> tmp_ver_err.txt || ver_ok=1
[[ $(wc -l < tmp_ver.txt) -eq 1 ]] && grep -q '^PLINK v1\.9' tmp_ver.txt && [[ ! -s tmp_ver_err.txt ]] || ver_ok=1
for v in -v -V; do
    plink $v 2>&1 | cmp -s - tmp_ver.txt || ver_ok=1
done
rm -f tmp_verx.* plink.log
plink --bfile tmp_cc --freq --version --out tmp_verx > tmp_ver2.txt 2>&1 || ver_ok=1
cmp -s tmp_ver.txt tmp_ver2.txt || ver_ok=1
! ls tmp_verx.* > /dev/null 2>&1 && [[ ! -e plink.log ]] || ver_ok=1
check "--version (and -v/-V) prints only the version line, ignoring other flags" $ver_ok

# --- --silent ---

# Nothing goes to stdout, but the outputs and log match a normal run, and
# warnings and errors still go to stderr.
si_ok=0
plink --bfile tmp_cc --freq --maf 0.1 --make-bed --seed 3 --out tmp_ns > tmp_ns.stdout 2> tmp_ns.stderr
plink --bfile tmp_cc --freq --maf 0.1 --make-bed --seed 3 --silent --out tmp_si > tmp_si.stdout 2> tmp_si.stderr || si_ok=1
[[ -s tmp_ns.stdout && ! -s tmp_si.stdout ]] || si_ok=1
cmp -s tmp_ns.stderr tmp_si.stderr || si_ok=1
# (--freq with --maf alone triggers an order-of-operations warning.)
plink --bfile tmp_cc --freq --maf 0.1 --silent --out tmp_siw > tmp_siw.stdout 2> tmp_siw.stderr || si_ok=1
[[ ! -s tmp_siw.stdout ]] && grep -q '^Warning' tmp_siw.stderr || si_ok=1
for e in frq bed bim fam; do
    test -s tmp_ns.$e && cmp -s tmp_ns.$e tmp_si.$e || si_ok=1
done
cmp -s <(log_body tmp_ns.log tmp_ns) <(log_body tmp_si.log tmp_si) || si_ok=1
cmp -s <(flags_plus tmp_ns.log tmp_ns --silent) <(flags_plus tmp_si.log tmp_si) || si_ok=1
plink --bfile tmp_nonexistent --freq --silent --out tmp_sie > tmp_sie.stdout 2> tmp_sie.stderr && si_ok=1
[[ ! -s tmp_sie.stdout ]] && grep -q '^Error' tmp_sie.stderr || si_ok=1
check "--silent suppresses only stdout" $si_ok

# --- --debug ---

# Only the logging method changes, so the outputs, console text and log are
# the same as without it.
db_ok=0
plink --bfile tmp_cc --freq --maf 0.1 --make-bed --seed 3 --debug --out tmp_db > tmp_db.stdout 2> tmp_db.stderr || db_ok=1
for e in frq bed bim fam; do
    test -s tmp_db.$e && cmp -s tmp_ns.$e tmp_db.$e || db_ok=1
done
cmp -s <(sed 's/tmp_ns/OUT/g' tmp_ns.stdout | grep -v '^  --') <(sed 's/tmp_db/OUT/g' tmp_db.stdout | grep -v '^  --') || db_ok=1
cmp -s <(log_body tmp_ns.log tmp_ns) <(log_body tmp_db.log tmp_db) || db_ok=1
cmp -s <(flags_plus tmp_ns.log tmp_ns --debug) <(flags_plus tmp_db.log tmp_db) || db_ok=1
# An error is logged too.
plink --bfile tmp_cc --maf 2 --make-bed --debug --out tmp_dbe > /dev/null 2>&1 && db_ok=1
grep -q '^Error' tmp_dbe.log || db_ok=1
check "--debug leaves outputs and log contents unchanged" $db_ok

# --- --script ---

# Flags in a script file, split across lines with tabs, extra spaces and a
# blank line, act as if they were on the command line in place of
# "--script <file>", including in the log.  A flag given both ways is an
# error rather than an override.
sc_ok=0
printf -- '--maf\t0.1\n\n   --geno 0.03   --hwe 0.001\n--freq\n' > tmp_script1.txt
printf -- '--bfile tmp_cc --maf 0.1\t--geno 0.03\n--hwe 0.001 --freq --make-bed --out tmp_scb\n' > tmp_script2.txt
plink --bfile tmp_cc --maf 0.1 --geno 0.03 --hwe 0.001 --freq --make-bed --out tmp_scd > /dev/null 2>&1
plink --bfile tmp_cc --script tmp_script1.txt --make-bed --out tmp_sca > /dev/null 2>&1 || sc_ok=1
plink --script tmp_script2.txt > /dev/null 2>&1 || sc_ok=1
for s in sca scb; do
    for e in frq bim fam; do
        test -s tmp_scd.$e && cmp -s tmp_scd.$e tmp_$s.$e || sc_ok=1
    done
    cmp -s <(log_flags tmp_scd.log | sed 's/tmp_scd/OUT/') <(log_flags tmp_$s.log | sed "s/tmp_$s/OUT/") || sc_ok=1
done
plink --bfile tmp_cc --script tmp_script1.txt --maf 0.2 --make-bed --out tmp_scx > tmp_scx.out 2>&1 && sc_ok=1
grep -q 'Duplicate --maf flag' tmp_scx.out || sc_ok=1
plink --script tmp_script1.txt --script tmp_script2.txt > tmp_scx.out 2>&1 && sc_ok=1
grep -q 'Multiple --script flags' tmp_scx.out || sc_ok=1
plink --script tmp_nonexistent.txt > tmp_scx.out 2>&1 && sc_ok=1
grep -q 'Failed to open tmp_nonexistent.txt' tmp_scx.out || sc_ok=1
check "--script expands a file of flags in place" $sc_ok

# --- --rerun ---

# The flags in a log are rerun; flags also on the command line (before or
# after --rerun) replace the logged ones, with all their parameters.
re_ok=0
plink --bfile tmp_cc --maf 0.1 --geno 0.03 --make-bed --out tmp_re1 > /dev/null 2>&1
plink --bfile tmp_cc --maf 0.25 --geno 0.03 --make-bed --out tmp_re0 > /dev/null 2>&1
cp tmp_re1.bim tmp_re1_orig.bim
plink --rerun tmp_re1.log > /dev/null 2>&1 || re_ok=1
cmp -s tmp_re1_orig.bim tmp_re1.bim || re_ok=1
plink --rerun tmp_re1.log --out tmp_re2 > /dev/null 2>&1 || re_ok=1
cmp -s tmp_re1.bim tmp_re2.bim && cmp -s tmp_re1.fam tmp_re2.fam || re_ok=1
cmp -s <(log_flags tmp_re1.log | sed 's/tmp_re1/OUT/') <(log_flags tmp_re2.log | sed 's/tmp_re2/OUT/') || re_ok=1
# The override must matter.
! cmp -s tmp_re0.bim tmp_re1.bim || re_ok=1
plink --rerun tmp_re1.log --maf 0.25 --out tmp_re3 > /dev/null 2>&1 || re_ok=1
plink --maf 0.25 --out tmp_re4 --rerun tmp_re1.log > /dev/null 2>&1 || re_ok=1
for s in re3 re4; do
    cmp -s tmp_re0.bim tmp_$s.bim || re_ok=1
    cmp -s <(log_flags tmp_re0.log | sed 's/tmp_re0/OUT/') <(log_flags tmp_$s.log | sed "s/tmp_$s/OUT/") || re_ok=1
done
# A log written from a --script run lists the expanded flags, so it reruns
# without the script.
plink --rerun tmp_sca.log --out tmp_re5 > /dev/null 2>&1 || re_ok=1
cmp -s tmp_scd.bim tmp_re5.bim && cmp -s tmp_scd.frq tmp_re5.frq || re_ok=1
check "--rerun repeats a logged command, with command-line overrides" $re_ok

# With no parameter, plink.log in the working directory is rerun.
rd_ok=0
mkdir -p tmp_rerun_dir
(
    cd tmp_rerun_dir || exit 1
    rm -f plink.*
    plink --bfile ../tmp_cc --maf 0.1 --geno 0.03 --make-bed > /dev/null 2>&1 || exit 1
    mv plink.bim first.bim
    plink --rerun > /dev/null 2>&1 || exit 1
    cmp -s first.bim plink.bim
) || rd_ok=1
check "--rerun defaults to plink.log" $rd_ok

# Older log formats: early PLINK 1.9 builds' "<n> arguments: ..." second line,
# and PLINK 1.07's tab-indented "Options in effect:" block.
ro_ok=0
{
    echo "PLINK v1.90b1"
    echo "9 arguments: --bfile tmp_cc --maf 0.1 --geno 0.03 --make-bed --out tmp_reo"
} > tmp_old_format.log
plink --rerun tmp_old_format.log > /dev/null 2>&1 || ro_ok=1
cmp -s tmp_re1.bim tmp_reo.bim || ro_ok=1
check "--rerun reads early PLINK 1.9 logs" $ro_ok
if [[ -x "$plink107" ]]; then
    r7_ok=0
    "$plink107" --noweb --bfile tmp_cc --maf 0.1 --geno 0.03 --make-bed --out tmp_re107 > /dev/null 2>&1
    grep -q "^	--maf 0.1" tmp_re107.log || r7_ok=1
    plink --rerun tmp_re107.log --out tmp_re7 > /dev/null 2>&1 || r7_ok=1
    cmp -s tmp_re1.bim tmp_re7.bim && cmp -s tmp_re1.fam tmp_re7.fam || r7_ok=1
    check "--rerun reads PLINK 1.07 logs" $r7_ok
else
    echo "note: PLINK 1.07 not found at $plink107; skipping its --rerun check"
fi

# Malformed or missing logs.
rx_ok=0
plink --rerun tmp_re1.log --rerun tmp_re1.log > tmp_rx.out 2>&1 && rx_ok=1
grep -q 'Duplicate --rerun flag' tmp_rx.out || rx_ok=1
plink --rerun tmp_nonexistent.log > tmp_rx.out 2>&1 && rx_ok=1
grep -q 'Failed to open --rerun file' tmp_rx.out || rx_ok=1
: > tmp_empty.log
plink --rerun tmp_empty.log > tmp_rx.out 2>&1 && rx_ok=1
grep -q 'Empty log file for --rerun' tmp_rx.out || rx_ok=1
echo "PLINK v1.9" > tmp_one_line.log
plink --rerun tmp_one_line.log > tmp_rx.out 2>&1 && rx_ok=1
grep -q 'Only one line in --rerun log file' tmp_rx.out || rx_ok=1
printf 'PLINK v1.9\nno flags here\n' > tmp_no_flags.log
plink --rerun tmp_no_flags.log > tmp_rx.out 2>&1 && rx_ok=1
grep -q 'Invalid log file for --rerun' tmp_rx.out || rx_ok=1
check "--rerun rejects duplicate flags and missing or malformed logs" $rx_ok

# --- --perm-batch-size ---

pb_ok=0
for v in 0 -1 1.5 abc; do
    plink --bfile tmp_cc --assoc mperm=10 --perm-batch-size $v --out tmp_pbx > tmp_pbx.out 2>&1 && pb_ok=1
    grep -q "Invalid --perm-batch-size parameter '$v'" tmp_pbx.out || pb_ok=1
done
plink --bfile tmp_cc --assoc mperm=10 --perm-batch-size --out tmp_pbx > /dev/null 2>&1 && pb_ok=1
check "--perm-batch-size rejects invalid parameters" $pb_ok

# With one thread, case/control permutations are generated one at a time
# from a single random stream, so batch boundaries (including partial final
# batches, batches of 1, and one batch larger than the permutation count)
# can't change the results.  This also holds for --linear and --qfam, which
# permute in the same order.
pi_ok=0
pi_cmds=(
    "--bfile tmp_cc --assoc mperm=300 --set-test --make-set tmp_sets.txt|assoc.set.mperm"
    "--bfile tmp_cc --model mperm=300 --set-test --make-set tmp_sets.txt|model.best.set.mperm"
    "--bfile tmp_cc --logistic mperm=200|assoc.logistic.mperm"
    "--bfile tmp_cc --logistic perm|assoc.logistic.perm"
    "--bfile tmp_qt --linear mperm=200|assoc.linear.mperm"
    "--bfile tmp_qt --linear mperm=200 --set-test --make-set tmp_sets.txt|assoc.linear.set.mperm"
    "--bfile tmp_fq --qfam mperm=300|qfam.within.perm"
)
for spec in "${pi_cmds[@]}"; do
    cmd=${spec%|*}
    ext=${spec#*|}
    rm -f tmp_pb_*.$ext
    for b in 512 1 7 64 5000; do
        plink $cmd --seed 5 --threads 1 --perm-batch-size $b --out tmp_pb_$b > /dev/null 2>&1
    done
    test -s tmp_pb_512.$ext || pi_ok=1
    for b in 1 7 64 5000; do
        cmp -s tmp_pb_512.$ext tmp_pb_$b.$ext || pi_ok=1
    done
done
check "--perm-batch-size doesn't change single-threaded case/control, --linear or --qfam permutation results" $pi_ok

# QT --assoc permutes in sample-major order within each batch, so the batch
# size changes which permutations are drawn.  The results must still be
# reproducible for a fixed batch size and thread count, and statistically
# equivalent across batch sizes.
pq_ok=0
for t in 1 4; do
    for b in 7 100; do
        plink --bfile tmp_qt --assoc mperm=2000 --seed 5 --threads $t --perm-batch-size $b --out tmp_pq_a > /dev/null 2>&1
        plink --bfile tmp_qt --assoc mperm=2000 --seed 5 --threads $t --perm-batch-size $b --out tmp_pq_b > /dev/null 2>&1
        test -s tmp_pq_a.qassoc.mperm && cmp -s tmp_pq_a.qassoc.mperm tmp_pq_b.qassoc.mperm || pq_ok=1
    done
done
plink --bfile tmp_qt --assoc mperm=2000 --seed 5 --threads 1 --perm-batch-size 7 --out tmp_pq_7 > /dev/null 2>&1
plink --bfile tmp_qt --assoc mperm=2000 --seed 5 --threads 1 --out tmp_pq_d > /dev/null 2>&1
! cmp -s tmp_pq_7.qassoc.mperm tmp_pq_d.qassoc.mperm || pq_ok=1
# Two independent EMP1 estimates from 2000 permutations differ with a
# standard deviation of at most 0.016; allow 5 of those for the worst of the
# 400 variants.
paste tmp_pq_7.qassoc.mperm tmp_pq_d.qassoc.mperm | awk 'NR > 1 {
    d = $3 - $7
    if (d < 0) d = -d
    sum += d
    if (d > max) max = d
    ++n
} END { exit !(n == 400 && sum / n < 0.015 && max < 0.08) }' || pq_ok=1
check "QT --assoc permutations are reproducible and batch-size-independent in distribution" $pq_ok

# QT --assoc --set-test.  Its worker threads used to read an unset thread
# count and skip every variant, so each permutation's set statistics came
# from uninitialized memory (usually zeroes, giving EMP1 = 1/(N+1) for every
# set with a significant variant).  On this null phenotype, EMP1 must be
# close to --linear's (the same t-statistics, permuted differently), at any
# batch size.
qs_ok=0
plink --bfile tmp_qt --linear mperm=2000 --set-test --make-set tmp_sets.txt --seed 5 --out tmp_qs_lin > /dev/null 2>&1
test -s tmp_qs_lin.assoc.linear.set.mperm || qs_ok=1
for b in 1 7 64 512 5000; do
    for t in 1 4; do
        rm -f tmp_qs.qassoc.set.mperm
        plink --bfile tmp_qt --assoc mperm=2000 --set-test --make-set tmp_sets.txt --seed 5 --threads $t --perm-batch-size $b --out tmp_qs > /dev/null 2>&1
        paste tmp_qs.qassoc.set.mperm tmp_qs_lin.assoc.linear.set.mperm | awk 'NR > 1 {
            ++n
            # SET and NSIG must match.
            if (($1 != $7) || ($3 != $9)) exit 1
            if ($3 != 0) {
                ++nsig
                d = $5 - $11
                if (d < 0) d = -d
                if (d > 0.06) exit 1
            }
        } END { exit !((n == 8) && (nsig >= 4)) }' || qs_ok=1
    done
done
check "QT --assoc --set-test permutation results match --linear's" $qs_ok

exit $failed

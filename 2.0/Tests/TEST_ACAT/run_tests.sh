#!/bin/bash

# "--acat-file" / "--set-list": the Cauchy combination of per-variant p-values,
# one p-value per set.
#
# Most of this runs against a hand-written results file rather than a --glm
# output, so that the expected values are exact identities of the test rather
# than whatever --glm happens to produce.

set -exo pipefail

BUILD=$1
EXTRA1=$2
EXTRA2=$3

# 1. A results file with p-values chosen to pin down each branch.
cat > tmp_res.txt << 'RESEOF'
#ID	TEST	P
v1	ADD	0.01
v2	ADD	0.01
v3	ADD	0.01
v4	ADD	0.01
v5	ADD	0.01
solo	ADD	0.037
mix1	ADD	0.9
mix2	ADD	0.004
mix3	ADD	0.25
tiny1	ADD	1e-300
tiny2	ADD	0.5
one1	ADD	1
one2	ADD	0.2
other	DOM	0.5
RESEOF

cat > tmp_sets.txt << 'SETEOF'
IDENTICAL 1 100 v1,v2,v3,v4,v5
SOLO 1 200 solo
MIXED 1 300 mix1,mix2,mix3
TINY 1 400 tiny1,tiny2
WITHONE 1 500 one1,one2
ABSENT 1 600 nosuchvariant
SETEOF

$BUILD/plink2 $EXTRA1 $EXTRA2 --acat-file tmp_res.txt test=ADD --set-list tmp_sets.txt --out tmp_a

# ABSENT has no member with a p-value, so it is dropped, not reported as 1.
if [[ $(grep -c '^ABSENT' tmp_a.acat) != 0 ]]; then
    echo "A set with no usable member was still reported."
    exit 1
fi
if [[ $(grep -cv '^#' tmp_a.acat) != 5 ]]; then
    echo "Expected 5 sets in the output."
    exit 1
fi

field() {
    awk -v s="$1" -v c="$2" '$1 == s {print $c}' tmp_a.acat
}
close_to() {
    awk -v got="$1" -v want="$2" -v tol="$3" -v lbl="$4" 'BEGIN {
        d = got - want; if (d < 0) d = -d;
        rel = d / (want < 0 ? -want : want);
        if (rel > tol) { print lbl ": got " got ", want " want; exit 1 }
    }'
}

# 2. The combination of k identical p-values is that p-value.  This is exact,
#    and independent of the weights, so it is the sharpest check available.
close_to "$(field IDENTICAL 6)" 0.01 1e-5 "IDENTICAL"

# 3. A one-member set returns that member's p-value.
close_to "$(field SOLO 6)" 0.037 1e-5 "SOLO"

# 4. NVAR counts what the set names; NVAR_TESTED counts what was found.
if [[ $(field MIXED 4) != 3 || $(field MIXED 5) != 3 ]]; then
    echo "MIXED variant counts wrong."
    exit 1
fi

# 5. Against the definition, computed independently here.  Equal weights, so
#    p = 0.5 - atan(mean of tan((0.5 - p_i) * pi)) / pi.
expected_mixed=$(awk 'BEGIN {
    pi = atan2(0, -1);
    n = split("0.9 0.004 0.25", p, " ");
    t = 0;
    for (i = 1; i <= n; ++i) { x = (0.5 - p[i]) * pi; t += sin(x) / cos(x) }
    r = t / n;
    if (r > 1) { print atan2(1.0 / r, 1) / pi } else { print 0.5 - atan2(r, 1) / pi }
}')
close_to "$(field MIXED 6)" "$expected_mixed" 1e-5 "MIXED"

# 6. A p-value small enough to underflow the naive formula must not destroy the
#    set: the combination is driven by it, and stays finite.
tiny=$(field TINY 6)
awk -v v="$tiny" 'BEGIN { if (!(v > 0) || v > 1e-299) { print "TINY: got " v; exit 1 } }'

# 7. A p-value of exactly 1 must not send the statistic to -infinity.
withone=$(field WITHONE 6)
awk -v v="$withone" 'BEGIN { if (!(v > 0) || !(v <= 1)) { print "WITHONE: got " v; exit 1 } }'

# 8. test= actually selects.  Dropping it leaves two rows for 'other' and the
#    DOM row would otherwise be pulled in; here we just check the flag is
#    required to be consistent.
$BUILD/plink2 $EXTRA1 $EXTRA2 --acat-file tmp_res.txt test=DOM --set-list tmp_sets.txt --out tmp_dom
if [[ $(grep -cv '^#' tmp_dom.acat) != 0 ]]; then
    echo "test=DOM should leave every set without a usable member."
    exit 1
fi

# 9. --set-list is required, and --acat-params only makes sense alongside.
fails() {
    if "$@" > /dev/null 2>&1; then
        echo "expected failure: $*"
        exit 1
    fi
}
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --acat-file tmp_res.txt --out tmp_bad
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --acat-params 1 25 --out tmp_bad
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --acat-file tmp_res.txt test=ADD --set-list tmp_sets.txt --acat-params 0 25 --out tmp_bad

# 10. End to end on a real --glm output, with the frequency-based weights.
printf '500 null 0.05 0.95 1 1\n' > tmp_sim.txt
plink --simulate tmp_sim.txt --make-bed --out tmp_g --seed 11 > /dev/null
{ echo "#FID IID QT"; awk '{print $1, $2, (NR % 89) / 89.0}' tmp_g.fam; } > tmp_qt.pheno
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_g --pheno tmp_qt.pheno --glm allow-no-covars cols=+a1freq --out tmp_gl
awk '{ i = int((NR - 1) / 10);
       ids[i] = (ids[i] ? ids[i] "," : "") $2; chr[i] = $1;
       if (!(i in pos)) pos[i] = $4 }
     END { for (j = 0; j <= i; ++j) print "S" j, chr[j], pos[j], ids[j] }' tmp_g.bim > tmp_gsets.txt
$BUILD/plink2 $EXTRA1 $EXTRA2 --acat-file tmp_gl.QT.glm.linear test=ADD --set-list tmp_gsets.txt --out tmp_ge
if [[ $(grep -cv '^#' tmp_ge.acat) != 50 ]]; then
    echo "Expected 50 sets from the --glm output."
    exit 1
fi
# Every reported p-value must be a probability.
awk '!/^#/ { if (!($6 > 0) || $6 > 1) { print "out-of-range p: " $0; exit 1 } }' tmp_ge.acat

echo "--acat-file tests passed."

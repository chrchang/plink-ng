# Recomputes the deterministic parts of a PLINK 1.9 case/control --assoc
# set test report: each set's size, significant variant count (NSIG),
# selected variant count (ISIG) and selected variants (SNPS), printed as
# "SET NSNP NSIG ISIG SNPS" (SNPS = NA when nothing is significant).  With
# ldset=1, prints the --set-r2 write (.ldset) contents instead: for each set
# with a significant variant, each member's in-set partners with r^2 >= r2max.
#
# Inputs, in order:
#   1. --make-set file (CHR START END NAME, inclusive; lines sharing a name
#      form one set);
#   2. the .bim file;
#   3. the --assoc report (CHISQ in column 8);
#   4. --r2 table output covering every same-chromosome pair.
# Variables: thr (the 1df chi-square value with p-value --set-p), lambda
# (--set-test-lambda; default 1), r2max (--set-r2; default 0.5), set_max
# (--set-max; default 5).
#
# A variant is significant when its chi-square divided by lambda is at least
# thr.  Significant variants are visited in decreasing chi-square order and
# selected unless their r^2 with an already-selected variant (same
# chromosome only) is at least r2max, until set_max are selected.
#
# Exits with status 2 if a chi-square or r^2 value is close enough to its
# cutoff, or two significant chi-squares close enough to each other, that
# rounding in the input files could change the result.

BEGIN {
    if (lambda == "") lambda = 1
    if (r2max == "") r2max = 0.5
    if (set_max == "") set_max = 5
}

FILENAME == ARGV[1] {
    if (!($4 in set_seen)) {
        set_seen[$4] = 1
        set_name[++set_ct] = $4
    }
    ++range_ct[$4]
    range_chr[$4, range_ct[$4]] = $1
    range_lo[$4, range_ct[$4]] = $2
    range_hi[$4, range_ct[$4]] = $3
    next
}

FILENAME == ARGV[2] {
    ++var_ct
    var_id[var_ct] = $2
    var_chr[var_ct] = $1
    var_pos[var_ct] = $4
    next
}

FILENAME == ARGV[3] {
    if (FNR > 1) {
        chisq[$2] = $8
        d = $8 / lambda - thr
        if ((d < 0.005 * thr) && (d > -0.005 * thr)) close_call = 1
    }
    next
}

FNR > 1 {
    r2[$3, $6] = $7
    r2[$6, $3] = $7
    d = $7 - r2max
    if ((d < 2e-6) && (d > -2e-6)) close_call = 1
}

function in_set(s, v,    k) {
    for (k = 1; k <= range_ct[s]; ++k) {
        if ((var_chr[v] == range_chr[s, k]) && (var_pos[v] >= range_lo[s, k]) && (var_pos[v] <= range_hi[s, k])) return 1
    }
    return 0
}

function conflict(a, b) {
    return ((a, b) in r2) && (r2[a, b] >= r2max)
}

END {
    for (si = 1; si <= set_ct; ++si) {
        s = set_name[si]
        member_ct = 0
        sig_ct = 0
        for (v = 1; v <= var_ct; ++v) {
            if (!in_set(s, v)) continue
            member[++member_ct] = var_id[v]
            if (chisq[var_id[v]] / lambda >= thr) sig[++sig_ct] = var_id[v]
        }
        if (ldset) {
            # (only sets with a significant variant are written)
            if (!sig_ct) continue
            for (a = 1; a <= member_ct; ++a) {
                line = ""
                for (b = 1; b <= member_ct; ++b) {
                    if ((a != b) && conflict(member[a], member[b])) line = line member[b] " "
                }
                if (line != "") print s, member[a], line
            }
            continue
        }
        # insertion sort, decreasing chi-square
        for (a = 2; a <= sig_ct; ++a) {
            x = sig[a]
            for (b = a - 1; (b >= 1) && (chisq[sig[b]] < chisq[x]); --b) sig[b + 1] = sig[b]
            sig[b + 1] = x
        }
        for (a = 2; a <= sig_ct; ++a) {
            if (chisq[sig[a]] >= 0.999 * chisq[sig[a - 1]]) close_call = 1
        }
        sel_ct = 0
        for (a = 1; (a <= sig_ct) && (sel_ct < set_max); ++a) {
            ok = 1
            for (b = 1; b <= sel_ct; ++b) {
                if (conflict(sel[b], sig[a])) {
                    ok = 0
                    break
                }
            }
            if (ok) sel[++sel_ct] = sig[a]
        }
        snps = sel_ct? sel[1] : "NA"
        for (b = 2; b <= sel_ct; ++b) snps = snps "|" sel[b]
        print s, member_ct, sig_ct, sel_ct, snps
    }
    if (close_call) exit 2
}

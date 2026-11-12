# Recounts transmitted (T) and untransmitted (U) A1 alleles from a .ped, and
# checks a PLINK 1.9 --tdt report's T, U, OR and CHISQ columns against them.
# A trio counts at a variant only when the father, mother and affected child
# are all genotyped there and the child's genotype is Mendelian-consistent.
#
# Usage: awk -f check_tdt.awk <.tdt> <.ped>
# Exits 1, printing the first few mismatches to stderr, on any mismatch.

function close_enough(a, b,    d, m) {
    if (a == "NA" || b == "NA") return a == b
    d = a - b
    if (d < 0) d = -d
    m = (a < 0)? -a : a
    if (b > m) m = b
    return d <= 1e-3 * m + 1e-6
}

function fail(msg) {
    if (++bad <= 5) print msg > "/dev/stderr"
}

# Number of A1 alleles in the genotype at variant j of the given .ped line, or
# -1 if missing.
function a1_ct(line_idx, j,    x, y) {
    x = geno[line_idx, 2 * j - 1]
    y = geno[line_idx, 2 * j]
    if (x == "0" || y == "0") return -1
    return (x == a1[j]) + (y == a1[j])
}

FNR == 1 { ++file_idx }

# .tdt: CHR SNP BP A1 A2 T U OR CHISQ P
file_idx == 1 {
    if (FNR == 1) next
    ++var_ct
    snp[var_ct] = $2
    a1[var_ct] = $4
    rep_t[var_ct] = $6
    rep_u[var_ct] = $7
    rep_or[var_ct] = $8
    rep_chisq[var_ct] = $9
    next
}

{
    ++line_ct
    iid[line_ct] = $2
    idx[$2] = line_ct
    pat[line_ct] = $3
    mat[line_ct] = $4
    pheno[line_ct] = $6
    for (i = 7; i <= NF; ++i) geno[line_ct, i - 6] = $i
}

END {
    for (j = 1; j <= var_ct; ++j) t[j] = u[j] = 0
    for (c = 1; c <= line_ct; ++c) {
        if (pheno[c] != 2 || !(pat[c] in idx) || !(mat[c] in idx)) continue
        f = idx[pat[c]]
        m = idx[mat[c]]
        for (j = 1; j <= var_ct; ++j) {
            gc = a1_ct(c, j)
            gf = a1_ct(f, j)
            gm = a1_ct(m, j)
            if (gc < 0 || gf < 0 || gm < 0) continue
            # Mendel check: the child's A1 count must be (0 or 1 from the
            # father) + (0 or 1 from the mother), as each genotype allows.
            lo = (gf == 2) + (gm == 2)
            hi = (gf > 0) + (gm > 0)
            if (gc < lo || gc > hi) continue
            if (gf == 1 && gm == 1) {
                t[j] += gc
                u[j] += 2 - gc
            } else if (gf == 1) {
                x = gc - (gm == 2)
                t[j] += x
                u[j] += 1 - x
            } else if (gm == 1) {
                x = gc - (gf == 2)
                t[j] += x
                u[j] += 1 - x
            }
        }
    }
    for (j = 1; j <= var_ct; ++j) {
        if (rep_t[j] != t[j] || rep_u[j] != u[j]) {
            fail(snp[j] ": T/U " rep_t[j] "/" rep_u[j] " vs " t[j] "/" u[j])
            continue
        }
        if (t[j] + u[j] == 0) continue
        if (!close_enough(rep_chisq[j], (t[j] - u[j]) ^ 2 / (t[j] + u[j]))) fail(snp[j] ": CHISQ " rep_chisq[j])
        if (u[j] > 0 && !close_enough(rep_or[j], t[j] / u[j])) fail(snp[j] ": OR " rep_or[j])
        ++checked
    }
    if (checked < var_ct / 2) fail("only " checked + 0 " of " var_ct " variants checked")
    exit (bad > 0)
}

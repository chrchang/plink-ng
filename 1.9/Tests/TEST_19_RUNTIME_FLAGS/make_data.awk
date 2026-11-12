# Writes a .ped of n samples over m A/C variants, with allele frequencies
# spread between 0.03 and 0.5 and about 2% missing calls.  Sexes alternate.
# Phenotypes are case/control (about half cases), or with qt=1, quantitative
# values in [0, 10).  With fam=1, samples come in families of four (two
# founders and two children, who inherit one random allele from each
# parent), so n should be a multiple of 4; otherwise they are unrelated.
#
# Usage: awk -v n=<n> -v m=<m> -v seed=<s> [-v qt=1] [-v fam=1]
#            -f make_data.awk > x.ped

function rnd() {
    s = (s * 69069 + 1) % 4294967296
    return s / 4294967296
}

function allele(j) {
    return (rnd() < freq[j])? "C" : "A"
}

function pheno() {
    return qt? sprintf("%.3f", rnd() * 10) : 1 + (rnd() < 0.5)
}

BEGIN {
    s = seed
    for (j = 1; j <= m; ++j) freq[j] = 0.03 + 0.47 * rnd()
    for (i = 1; i <= n; ++i) {
        r = (i - 1) % 4
        if (!fam) {
            fid = "f" i
            pat = mat = 0
        } else {
            fid = "f" (1 + int((i - 1) / 4))
            if (r < 2) {
                pat = mat = 0
            } else {
                pat = "i" (i - r)
                mat = "i" (i - r + 1)
            }
        }
        line = fid " i" i " " pat " " mat " " (1 + (i % 2)) " " pheno()
        for (j = 1; j <= m; ++j) {
            if (fam && (r >= 2)) {
                # One allele from each parent (blanked below if missing).
                a1 = g[i - r, j, 1 + (rnd() < 0.5)]
                a2 = g[i - r + 1, j, 1 + (rnd() < 0.5)]
            } else {
                a1 = allele(j)
                a2 = allele(j)
            }
            g[i, j, 1] = a1
            g[i, j, 2] = a2
            if (rnd() < 0.02) {
                line = line " 0 0"
            } else {
                line = line " " a1 " " a2
            }
        }
        print line
    }
}

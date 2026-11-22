# Writes a .ped of n_fam families (two founders and two children, who inherit
# one random allele from each parent) followed by n_unrel unrelated samples,
# over n_var A/C variants with allele frequencies between 0.1 and 0.9.  Every
# fourth sample (by ID number) is affected; with qt=1, phenotypes are instead
# quantitative values.
#
# Usage: awk -v n_fam=<n> -v n_unrel=<n> -v n_var=<m> -v seed=<s> [-v qt=1]
#            -f make_families.awk > x.ped

function rnd() {
    s = (s * 69069 + 1) % 4294967296
    return s / 4294967296
}

function allele(j) {
    return (rnd() < freq[j])? "C" : "A"
}

function pheno(i) {
    return qt? sprintf("%.3f", rnd() * 10) : 1 + (i % 4 == 0)
}

BEGIN {
    s = seed
    for (j = 1; j <= n_var; ++j) freq[j] = 0.1 + 0.8 * rnd()
    id = 0
    for (f = 1; f <= n_fam; ++f) {
        for (p = 1; p <= 2; ++p) {
            ++id
            line[p] = "f" f " i" id " 0 0 " p " " pheno(id)
            for (j = 1; j <= n_var; ++j) {
                h[p, j, 1] = allele(j)
                h[p, j, 2] = allele(j)
                line[p] = line[p] " " h[p, j, 1] " " h[p, j, 2]
            }
            pid[p] = "i" id
        }
        print line[1]
        print line[2]
        for (c = 1; c <= 2; ++c) {
            ++id
            out = "f" f " i" id " " pid[1] " " pid[2] " " (1 + id % 2) " " pheno(id)
            for (j = 1; j <= n_var; ++j) out = out " " h[1, j, 1 + (rnd() < 0.5)] " " h[2, j, 1 + (rnd() < 0.5)]
            print out
        }
    }
    for (u = 1; u <= n_unrel; ++u) {
        ++id
        out = "u" u " i" id " 0 0 " (1 + id % 2) " " pheno(id)
        for (j = 1; j <= n_var; ++j) out = out " " allele(j) " " allele(j)
        print out
    }
}

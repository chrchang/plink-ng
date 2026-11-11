# Writes a .ped of n_trio trios (father, mother, affected child) plus one
# unrelated sample, over n_var biallelic A/C variants.  Children inherit one
# random allele from each parent.  About 2% of calls are missing, and about
# 1% of children's calls are replaced with a random genotype, which creates
# some Mendel errors.
#
# Usage: awk -v n_trio=<n> -v n_var=<m> -v seed=<s> -f make_trios.awk > x.ped

function rnd() {
    s = (s * 69069 + 1) % 4294967296
    return s / 4294967296
}

function allele(p) {
    return (rnd() < p)? "A" : "C"
}

function call(a, b) {
    return (rnd() < 0.02)? "0 0" : a " " b
}

BEGIN {
    s = seed
    for (j = 1; j <= n_var; ++j) freq[j] = 0.1 + 0.8 * rnd()
    for (t = 1; t <= n_trio; ++t) {
        fa = "fa" t
        mo = "mo" t
        ch = "ch" t
        fline = "f" t " " fa " 0 0 1 1"
        mline = "f" t " " mo " 0 0 2 1"
        cline = "f" t " " ch " " fa " " mo " " (1 + (t % 2)) " 2"
        for (j = 1; j <= n_var; ++j) {
            f1 = allele(freq[j]); f2 = allele(freq[j])
            m1 = allele(freq[j]); m2 = allele(freq[j])
            c1 = (rnd() < 0.5)? f1 : f2
            c2 = (rnd() < 0.5)? m1 : m2
            if (rnd() < 0.01) {
                c1 = allele(0.5)
                c2 = allele(0.5)
            }
            fline = fline " " call(f1, f2)
            mline = mline " " call(m1, m2)
            cline = cline " " call(c1, c2)
        }
        print fline
        print mline
        print cline
    }
    line = "u1 u1 0 0 1 1"
    for (j = 1; j <= n_var; ++j) line = line " " call(allele(freq[j]), allele(freq[j]))
    print line
}

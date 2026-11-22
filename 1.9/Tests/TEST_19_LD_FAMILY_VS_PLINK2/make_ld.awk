# Writes a .ped of n unrelated samples over m biallelic A/C variants with
# linkage disequilibrium: each block of 10 variants follows one of 4
# haplotypes per chromosome copy (2 bits).  Within a block, variants j with
# j % 3 == 2 carry bit 0, j % 3 == 0 carry bit 1, and j % 3 == 1 are 1 only
# for haplotype 3 (so they have r^2 of about 1/3 with each bit).  Alleles are
# flipped with 5% probability, and 1% of genotypes are missing.  Sex
# alternates; every third sample is a case.
#
# Usage: awk -v n=<samples> -v m=<variants> -v seed=<s> -f make_ld.awk > x.ped

function rnd() {
    s = (s * 69069 + 1) % 4294967296
    return int(s / 65536)
}

BEGIN {
    s = seed
    for (i = 1; i <= n; ++i) {
        line = "f" i " i" i " 0 0 " (1 + i % 2) " " (1 + (i % 3 == 0))
        for (j = 1; j <= m; ++j) {
            if (j % 10 == 1) {
                for (h = 1; h <= 2; ++h) hap[h] = rnd() % 4
            }
            if (rnd() % 100 < 1) {
                line = line " 0 0"
                continue
            }
            for (h = 1; h <= 2; ++h) {
                a = (j % 3 == 1)? (hap[h] == 3) : int(hap[h] / (1 + (j % 3 == 0))) % 2
                if (rnd() % 100 < 5) a = 1 - a
                line = line " " (a? "C" : "A")
            }
        }
        print line
    }
}

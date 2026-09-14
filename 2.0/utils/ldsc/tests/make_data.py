#!/usr/bin/env python3
"""Simulate LD Scores and GWAS summary statistics under the LDSC model.

For each variant j, with M variants and LD Score l2_j,

    var(z1_j) = intercept1 + N1 * h2_1 * l2_j / M
    var(z2_j) = intercept2 + N2 * h2_2 * l2_j / M
    cov(z1_j, z2_j) = gencov_intercept + sqrt(N1 N2) * rho_g * l2_j / M

which is the model the regression inverts, so the estimates land near the
parameters set here.  N/M is kept realistic (mean chi^2 of about 2) rather
than scaled up with M: the whole point of the two-step estimator is what it
does to the high-chi^2 tail, and an unrealistic mean chi^2 would put most of
the data in that tail.  Deliberately numpy-free: the test runners do not have
numpy.
"""
import math
import random
import sys

M = 6000
N1 = 3000.0
N2 = 2400.0
H2_1 = 0.30
H2_2 = 0.18
RHO_G = 0.10
INTERCEPT1 = 1.08
INTERCEPT2 = 1.04
GENCOV_INTERCEPT = 0.03
ALLELE_PAIRS = (('A', 'G'), ('C', 'T'), ('A', 'C'), ('G', 'T'))
# Partitioned fixture: three non-overlapping annotations with different
# per-variant heritability.
PART_M = (3000.0, 2000.0, 1000.0)
PART_H2 = (0.18, 0.09, 0.04)
PART_INTERCEPT = 1.05
# Overlapping fixture: a baseline annotation covering everything, plus two
# that overlap each other.  Written with .annot and .frq files, so the overlap
# correction has something to correct with.
OV_FRAC = (1.0, 0.40, 0.25)
OV_H2 = (0.15, 0.08, 0.05)


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    rand = random.Random(seed)
    ids = ['rs%d' % (j + 1) for j in range(M)]
    # Gamma-ish LD Scores: a sum of exponentials, so no numpy needed.
    l2 = [0.5 - 3.0 * (math.log(rand.random()) + math.log(rand.random()))
          for _ in range(M)]
    # The regression weight LD Scores are a noisy version of the same thing,
    # as they are in practice (sum of r^2 over a different variant set).
    w_ld = [l2[j] * rand.uniform(0.8, 1.2) for j in range(M)]

    z1 = []
    z2 = []
    for j in range(M):
        var1 = INTERCEPT1 + N1 * H2_1 * l2[j] / M
        var2 = INTERCEPT2 + N2 * H2_2 * l2[j] / M
        cov = GENCOV_INTERCEPT + math.sqrt(N1 * N2) * RHO_G * l2[j] / M
        # Cholesky of a 2x2, applied to two standard normals.
        a = math.sqrt(var1)
        b = cov / a
        c_sq = var2 - b * b
        c = math.sqrt(c_sq) if c_sq > 0.0 else 0.0
        u1 = rand.gauss(0.0, 1.0)
        u2 = rand.gauss(0.0, 1.0)
        z1.append(a * u1)
        z2.append(b * u1 + c * u2)

    with open('ldscores.ldscore', 'w') as f:
        f.write('#CHROM\tPOS\tID\tL2\n')
        for j in range(M):
            f.write('1\t%d\t%s\t%.6f\n' % (j + 1, ids[j], l2[j]))
    with open('w_ld.ldscore', 'w') as f:
        f.write('#CHROM\tPOS\tID\tL2\n')
        for j in range(M):
            f.write('1\t%d\t%s\t%.6f\n' % (j + 1, ids[j], w_ld[j]))
    with open('ldscores.ldscore.l2.M_5_50', 'w') as f:
        f.write('%d\n' % M)

    # One allele pair per variant, shared by both traits, as munged summary
    # statistics have.  (Test 6 is what exercises disagreeing orientations.)
    alleles = [ALLELE_PAIRS[rand.randrange(4)] for _ in range(M)]
    for name, z, n in (('trait1', z1, N1), ('trait2', z2, N2)):
        with open(name + '.sumstats', 'w') as f:
            f.write('SNP\tA1\tA2\tZ\tN\n')
            for j in range(M):
                a1, a2 = alleles[j]
                f.write('%s\t%s\t%s\t%.6f\t%g\n' % (ids[j], a1, a2, z[j], n))

    # trait2 again, with the effect alleles swapped on half the variants (Z
    # negated to match) and the strand flipped on a quarter.  Estimates must
    # not change.
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    with open('trait2_flipped.sumstats', 'w') as f:
        f.write('SNP\tA1\tA2\tZ\tN\n')
        with open('trait2.sumstats') as src:
            src.readline()
            for j, line in enumerate(src):
                snp, a1, a2, z, n = line.split()
                if j % 2 == 0:
                    a1, a2 = a2, a1
                    z = '%.6f' % -float(z)
                if j % 4 == 1:
                    a1, a2 = complement[a1], complement[a2]
                f.write('%s\t%s\t%s\t%s\t%s\n' % (snp, a1, a2, z, n))

    # The same LD Scores split across three chromosomes, in ldsc's own naming.
    per_chr = (M + 2) // 3
    for chrom in (1, 2, 3):
        lo = (chrom - 1) * per_chr
        hi = min(chrom * per_chr, M)
        with open('chr_ref.%d.l2.ldscore' % chrom, 'w') as f:
            f.write('CHR\tSNP\tL2\n')
            for j in range(lo, hi):
                f.write('%d\t%s\t%.6f\n' % (chrom, ids[j], l2[j]))
        with open('chr_ref.%d.l2.M_5_50' % chrom, 'w') as f:
            f.write('%d\n' % (hi - lo))
        with open('chr_w.%d.l2.ldscore' % chrom, 'w') as f:
            f.write('CHR\tSNP\tL2\n')
            for j in range(lo, hi):
                f.write('%d\t%s\t%.6f\n' % (chrom, ids[j], w_ld[j]))

    # A partitioned LD Score fileset: one L2 column per annotation, with the
    # per-annotation variant counts in the .l2.M_5_50 file, plus a trait
    # simulated from tau_c = h2_c / M_c.
    part_l2 = []
    for j in range(M):
        part_l2.append([0.2 - 2.0 * (math.log(rand.random()) +
                                     math.log(rand.random()))
                        for _ in range(len(PART_M))])
    tau = [PART_H2[c] / PART_M[c] for c in range(len(PART_M))]
    part_z = []
    for j in range(M):
        var = PART_INTERCEPT + N1 * sum(part_l2[j][c] * tau[c]
                                        for c in range(len(PART_M)))
        part_z.append(math.sqrt(var) * rand.gauss(0.0, 1.0))
    with open('part_ref.l2.ldscore', 'w') as f:
        f.write('CHR\tSNP\tBP\tannotAL2\tannotBL2\tannotCL2\n')
        for j in range(M):
            f.write('1\t%s\t%d\t%.6f\t%.6f\t%.6f\n'
                    % (ids[j], j + 1, part_l2[j][0], part_l2[j][1],
                       part_l2[j][2]))
    with open('part_ref.l2.M_5_50', 'w') as f:
        f.write('%g %g %g\n' % PART_M)
    with open('part_w.l2.ldscore', 'w') as f:
        f.write('CHR\tSNP\tBP\tL2\n')
        for j in range(M):
            f.write('1\t%s\t%d\t%.6f\n'
                    % (ids[j], j + 1, sum(part_l2[j]) * rand.uniform(0.85,
                                                                     1.15)))
    with open('part_trait.sumstats', 'w') as f:
        f.write('SNP\tA1\tA2\tZ\tN\n')
        for j in range(M):
            a1, a2 = alleles[j]
            f.write('%s\t%s\t%s\t%.6f\t%g\n'
                    % (ids[j], a1, a2, part_z[j], N1))

    # Overlapping annotations, with the .annot and .frq files --overlap-annot
    # needs.  Only the common variants (5% < MAF < 50%) count, as in the
    # .l2.M_5_50 convention.
    ov_frq = [rand.uniform(0.01, 0.5) for _ in range(M)]
    ov_annot = []
    for j in range(M):
        ov_annot.append([1.0] + [1.0 if rand.random() < OV_FRAC[c] else 0.0
                                 for c in range(1, len(OV_FRAC))])
    common = [j for j in range(M) if 0.05 < ov_frq[j] < 0.95]
    ov_m = [sum(ov_annot[j][c] for j in common) for c in range(len(OV_FRAC))]
    ov_tau = [OV_H2[c] / ov_m[c] for c in range(len(OV_FRAC))]
    ov_l2 = []
    for j in range(M):
        row = []
        for c in range(len(OV_FRAC)):
            base = 0.1 - 2.0 * (math.log(rand.random()) +
                                math.log(rand.random()))
            row.append(base * (0.3 + 0.7 * ov_annot[j][c]))
        ov_l2.append(row)
    ov_z = []
    for j in range(M):
        var = 1.03 + N1 * sum(ov_l2[j][c] * ov_tau[c]
                              for c in range(len(OV_FRAC)))
        ov_z.append(math.sqrt(var) * rand.gauss(0.0, 1.0))
    with open('ov_ref.l2.ldscore', 'w') as f:
        f.write('CHR\tSNP\tBP\tbaseL2\tannotBL2\tannotCL2\n')
        for j in range(M):
            f.write('1\t%s\t%d\t%.6f\t%.6f\t%.6f\n'
                    % (ids[j], j + 1, ov_l2[j][0], ov_l2[j][1], ov_l2[j][2]))
    with open('ov_ref.l2.M_5_50', 'w') as f:
        f.write('%g %g %g\n' % tuple(ov_m))
    with open('ov_ref.annot', 'w') as f:
        f.write('CHR\tBP\tSNP\tCM\tbase\tannotB\tannotC\n')
        for j in range(M):
            f.write('1\t%d\t%s\t0\t%g\t%g\t%g\n'
                    % (j + 1, ids[j], ov_annot[j][0], ov_annot[j][1],
                       ov_annot[j][2]))
    with open('ov_ref.frq', 'w') as f:
        f.write('CHR\tSNP\tA1\tA2\tMAF\tNCHROBS\n')
        for j in range(M):
            f.write('1\t%s\tA\tG\t%.6f\t1000\n' % (ids[j], ov_frq[j]))
    with open('ov_w.l2.ldscore', 'w') as f:
        f.write('CHR\tSNP\tBP\tL2\n')
        for j in range(M):
            f.write('1\t%s\t%d\t%.6f\n'
                    % (ids[j], j + 1,
                       sum(ov_l2[j]) * rand.uniform(0.85, 1.15)))
    with open('ov_trait.sumstats', 'w') as f:
        f.write('SNP\tA1\tA2\tZ\tN\n')
        for j in range(M):
            a1, a2 = alleles[j]
            f.write('%s\t%s\t%s\t%.6f\t%g\n'
                    % (ids[j], a1, a2, ov_z[j], N1))

    print('%d variants' % M)


if __name__ == '__main__':
    main()

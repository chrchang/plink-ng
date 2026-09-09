#!/usr/bin/env python3
"""Write a trio dataset for --tucc, plus the pseudo-case/control answer.

Parental genotypes are drawn from an allele frequency, the child's alleles are
drawn from the parents, and then a fraction of the calls is overwritten with a
Mendel error or a missing call so that both of those paths are exercised.
"""
import random
import sys

FAMILY_CT = 40
AUTO_VARIANT_CT = 300


def geno_str(alt_ct):
    if alt_ct is None:
        return './.'
    return {0: '0/0', 1: '0/1', 2: '1/1'}[alt_ct]


def main():
    seed = int(sys.argv[1])
    rand = random.Random(seed)
    sample_ids = []
    for family_idx in range(FAMILY_CT):
        sample_ids.extend(['F%d_DAD' % family_idx, 'F%d_MOM' % family_idx,
                           'F%d_KID' % family_idx])
    # Two of the families are broken on purpose: one child has a single parent
    # in the dataset, and one is a founder.  Neither may show up in the output.
    with open('trios.psam', 'w') as f:
        f.write('#FID\tIID\tPAT\tMAT\tSEX\tPHENO1\n')
        for family_idx in range(FAMILY_CT):
            dad, mom, kid = ['F%d_%s' % (family_idx, s)
                             for s in ('DAD', 'MOM', 'KID')]
            fid = 'FAM%d' % family_idx
            f.write('%s\t%s\t0\t0\t1\t1\n' % (fid, dad))
            f.write('%s\t%s\t0\t0\t2\t1\n' % (fid, mom))
            kid_sex = 1 + (family_idx % 2)
            if family_idx == FAMILY_CT - 1:
                kid_pat, kid_mat = '0', '0'
            elif family_idx == FAMILY_CT - 2:
                kid_pat, kid_mat = dad, '0'
            else:
                kid_pat, kid_mat = dad, mom
            f.write('%s\t%s\t%s\t%s\t%d\t2\n' % (fid, kid, kid_pat, kid_mat,
                                                 kid_sex))

    rows = []
    # 'chr1' carries the variants --tucc must keep; chrX and chrMT must be
    # dropped, and so must the multiallelic variant.
    for variant_idx in range(AUTO_VARIANT_CT):
        rows.append(('1', variant_idx + 1, 'v%d' % variant_idx, 'A', 'G'))
    rows.append(('1', AUTO_VARIANT_CT + 1, 'vmulti', 'A', 'G,T'))
    rows.append(('X', 60000000, 'vx', 'A', 'G'))
    rows.append(('MT', 1, 'vmt', 'A', 'G'))

    lines = []
    expected = {}
    for chrom, pos, vid, ref, alt in rows:
        freq = rand.uniform(0.1, 0.9)
        calls = {}
        for family_idx in range(FAMILY_CT):
            dad_alleles = [int(rand.random() < freq) for _ in range(2)]
            mom_alleles = [int(rand.random() < freq) for _ in range(2)]
            transmitted = [rand.choice(dad_alleles), rand.choice(mom_alleles)]
            dad_geno = sum(dad_alleles)
            mom_geno = sum(mom_alleles)
            kid_geno = sum(transmitted)
            roll = rand.random()
            if roll < 0.05:
                kid_geno = None
            elif roll < 0.10:
                dad_geno = None
            elif roll < 0.15:
                # A Mendel error, when one is reachable.
                for candidate in (0, 1, 2):
                    lo = (dad_geno == 2) + (mom_geno == 2)
                    hi = (dad_geno != 0) + (mom_geno != 0)
                    if candidate < lo or candidate > hi:
                        kid_geno = candidate
                        break
            calls[family_idx] = (dad_geno, mom_geno, kid_geno)
        if alt == 'G' and chrom == '1':
            expected[vid] = calls
        geno_fields = []
        for family_idx in range(FAMILY_CT):
            dad_geno, mom_geno, kid_geno = calls[family_idx]
            geno_fields.extend([geno_str(dad_geno), geno_str(mom_geno),
                                geno_str(kid_geno)])
        lines.append('%s\t%d\t%s\t%s\t%s\t.\t.\t.\tGT\t%s' %
                     (chrom, pos, vid, ref, alt, '\t'.join(geno_fields)))

    with open('trios.vcf', 'w') as f:
        f.write('##fileformat=VCFv4.2\n')
        for chrom in ('1', 'X', 'MT'):
            f.write('##contig=<ID=%s>\n' % chrom)
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' %
                '\t'.join(sample_ids))
        for line in lines:
            f.write(line + '\n')

    # The answer, as ALT counts: for a complete Mendel-consistent trio the
    # pseudo-case is the child and the pseudo-control is the two alleles the
    # parents kept, which is (dad + mom) - child.
    with open('expected.txt', 'w') as f:
        for variant_idx in range(AUTO_VARIANT_CT):
            vid = 'v%d' % variant_idx
            calls = expected[vid]
            out = []
            for family_idx in range(FAMILY_CT - 2):
                dad_geno, mom_geno, kid_geno = calls[family_idx]
                case_geno, control_geno = 'NA', 'NA'
                if None not in (dad_geno, mom_geno, kid_geno):
                    lo = (dad_geno == 2) + (mom_geno == 2)
                    hi = (dad_geno != 0) + (mom_geno != 0)
                    if lo <= kid_geno <= hi:
                        case_geno = str(kid_geno)
                        control_geno = str(dad_geno + mom_geno - kid_geno)
                out.extend([case_geno, control_geno])
            f.write('%s\t%s\n' % (vid, '\t'.join(out)))


if __name__ == '__main__':
    main()

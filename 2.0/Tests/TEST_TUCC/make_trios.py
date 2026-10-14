#!/usr/bin/env python3
"""Write a trio dataset for --tucc, plus the pseudo-case/control answer.

Parental genotypes are drawn from an allele frequency, the child's alleles are
drawn from the parents, and then a fraction of the calls is overwritten with a
Mendel error or a missing call so that both of those paths are exercised.

Biallelic variants are checked as ALT counts; the multiallelic ones are
checked as allele pairs, since that is the only way the untransmitted alleles
are pinned down.
"""
import random
import sys

FAMILY_CT = 40
AUTO_VARIANT_CT = 300
MULTI_VARIANT_CT = 60
MULTI_ALLELE_CT = 3


def geno_str(alt_ct):
    if alt_ct is None:
        return './.'
    return {0: '0/0', 1: '0/1', 2: '1/1'}[alt_ct]


def pair_str(pair):
    if pair is None:
        return './.'
    return '%d/%d' % (min(pair), max(pair))


def transmissible(parent, allele):
    return allele in parent


def consistent(dad, mom, kid):
    """True iff one of the child's alleles can come from each parent."""
    return ((transmissible(dad, kid[0]) and transmissible(mom, kid[1])) or
            (transmissible(dad, kid[1]) and transmissible(mom, kid[0])))


def untransmitted(dad, mom, kid):
    """The pair of parental alleles the child did not receive."""
    def other(parent, transmitted):
        return parent[1] if parent[0] == transmitted else parent[0]
    if transmissible(dad, kid[0]) and transmissible(mom, kid[1]):
        return (other(dad, kid[0]), other(mom, kid[1]))
    return (other(dad, kid[1]), other(mom, kid[0]))


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
    for variant_idx in range(MULTI_VARIANT_CT):
        rows.append(('1', AUTO_VARIANT_CT + variant_idx + 1,
                     'vm%d' % variant_idx, 'A', 'G,T'))
    rows.append(('X', 60000000, 'vx', 'A', 'G'))
    rows.append(('MT', 1, 'vmt', 'A', 'G'))

    lines = []
    expected = {}
    expected_multi = {}
    for chrom, pos, vid, ref, alt in rows:
        is_multi = ',' in alt
        calls = {}
        if is_multi:
            weights = [rand.uniform(0.1, 0.9) for _ in range(MULTI_ALLELE_CT)]
            total = sum(weights)
            weights = [w / total for w in weights]

            def draw():
                roll = rand.random()
                acc = 0.0
                for allele_idx, weight in enumerate(weights):
                    acc += weight
                    if roll < acc:
                        return allele_idx
                return MULTI_ALLELE_CT - 1

            for family_idx in range(FAMILY_CT):
                dad = (draw(), draw())
                mom = (draw(), draw())
                kid = (rand.choice(dad), rand.choice(mom))
                roll = rand.random()
                if roll < 0.05:
                    kid = None
                elif roll < 0.10:
                    dad = None
                elif roll < 0.20:
                    # A Mendel error, when one is reachable.
                    for cand0 in range(MULTI_ALLELE_CT):
                        for cand1 in range(MULTI_ALLELE_CT):
                            if not consistent(dad, mom, (cand0, cand1)):
                                kid = (cand0, cand1)
                                break
                        else:
                            continue
                        break
                calls[family_idx] = (dad, mom, kid)
            expected_multi[vid] = calls
            geno_fields = []
            for family_idx in range(FAMILY_CT):
                dad, mom, kid = calls[family_idx]
                geno_fields.extend([pair_str(dad), pair_str(mom),
                                    pair_str(kid)])
        else:
            freq = rand.uniform(0.1, 0.9)
            for family_idx in range(FAMILY_CT):
                dad_alleles = [int(rand.random() < freq) for _ in range(2)]
                mom_alleles = [int(rand.random() < freq) for _ in range(2)]
                transmitted = [rand.choice(dad_alleles),
                               rand.choice(mom_alleles)]
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
            if chrom == '1':
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

    # The multiallelic answer, as allele-index pairs: the pseudo-case is the
    # child and the pseudo-control is each parent's non-transmitted allele.
    with open('expected_multi.txt', 'w') as f:
        for variant_idx in range(MULTI_VARIANT_CT):
            vid = 'vm%d' % variant_idx
            calls = expected_multi[vid]
            out = []
            for family_idx in range(FAMILY_CT - 2):
                dad, mom, kid = calls[family_idx]
                case_geno, control_geno = './.', './.'
                if None not in (dad, mom, kid) and consistent(dad, mom, kid):
                    case_geno = pair_str(kid)
                    control_geno = pair_str(untransmitted(dad, mom, kid))
                out.extend([case_geno, control_geno])
            f.write('%s\t%s\n' % (vid, '\t'.join(out)))


if __name__ == '__main__':
    main()

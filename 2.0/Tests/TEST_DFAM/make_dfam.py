#!/usr/bin/env python3
"""Write a --dfam fixture exercising all four kinds of group.

  * nuclear families whose genotyped children are all cases,
  * nuclear families with both case and control children,
  * sibships, i.e. samples sharing FID/PAT/MAT whose parents are not in the
    dataset,
  * unrelated founders.

Parental genotypes are drawn from a per-variant allele frequency and children's
alleles are drawn from the parents, with Mendel errors and missing calls mixed
in.
"""
import random
import sys

ALL_CASE_FAMILY_CT = 25
MIXED_FAMILY_CT = 25
SIBSHIP_CT = 20
UNRELATED_CT = 60
VARIANT_CT = 400


def geno_str(alt_ct):
    if alt_ct is None:
        return './.'
    return {0: '0/0', 1: '0/1', 2: '1/1'}[alt_ct]


def draw_child(rand, dad_geno, mom_geno):
    dad_alleles = [1, 1] if dad_geno == 2 else ([0, 1] if dad_geno == 1 else [0, 0])
    mom_alleles = [1, 1] if mom_geno == 2 else ([0, 1] if mom_geno == 1 else [0, 0])
    return rand.choice(dad_alleles) + rand.choice(mom_alleles)


def main():
    rand = random.Random(int(sys.argv[1]))
    samples = []  # (fid, iid, pat, mat, sex, pheno)
    # kind: ('family', dad_iid, mom_iid) for children, None for everyone else
    parents_of = {}

    def add(fid, iid, pat, mat, sex, pheno):
        samples.append((fid, iid, pat, mat, sex, pheno))

    for family_idx in range(ALL_CASE_FAMILY_CT + MIXED_FAMILY_CT):
        all_case = family_idx < ALL_CASE_FAMILY_CT
        fid = 'F%d' % family_idx
        dad, mom = '%s_D' % fid, '%s_M' % fid
        add(fid, dad, '0', '0', 1, 1)
        add(fid, mom, '0', '0', 2, 1)
        child_ct = 2 + (family_idx % 3)
        for child_idx in range(child_ct):
            iid = '%s_C%d' % (fid, child_idx)
            if all_case:
                pheno = 2
            else:
                pheno = 2 if child_idx == 0 else (1 if child_idx == 1 else
                                                  rand.choice([1, 2]))
            add(fid, iid, dad, mom, 1 + (child_idx % 2), pheno)
            parents_of[iid] = (dad, mom)

    for sibship_idx in range(SIBSHIP_CT):
        fid = 'S%d' % sibship_idx
        # The parents are named but absent from the dataset, which is what
        # makes this a sibship rather than a nuclear family.
        dad, mom = '%s_ABSENT_D' % fid, '%s_ABSENT_M' % fid
        sib_ct = 2 + (sibship_idx % 3)
        for sib_idx in range(sib_ct):
            pheno = 2 if sib_idx == 0 else (1 if sib_idx == 1 else
                                            rand.choice([1, 2]))
            add(fid, '%s_S%d' % (fid, sib_idx), dad, mom, 1 + (sib_idx % 2),
                pheno)
            parents_of['%s_S%d' % (fid, sib_idx)] = (dad, mom)

    for unrelated_idx in range(UNRELATED_CT):
        add('U%d' % unrelated_idx, 'U%d_I' % unrelated_idx, '0', '0',
            1 + (unrelated_idx % 2), 1 + (unrelated_idx % 3 == 0))

    with open('dfam.psam', 'w') as f:
        f.write('#FID\tIID\tPAT\tMAT\tSEX\tPHENO1\n')
        for row in samples:
            f.write('%s\t%s\t%s\t%s\t%d\t%d\n' % row)

    sample_ids = [row[1] for row in samples]
    id_set = set(sample_ids)
    lines = []
    for variant_idx in range(VARIANT_CT):
        freq = rand.uniform(0.15, 0.85)
        calls = {}
        # Founders (and absent parents' stand-ins) first, then children.
        for fid, iid, pat, mat, _sex, _pheno in samples:
            if iid in parents_of:
                continue
            calls[iid] = sum(1 for _ in range(2) if rand.random() < freq)
        for fid, iid, pat, mat, _sex, _pheno in samples:
            if iid not in parents_of:
                continue
            dad, mom = parents_of[iid]
            if dad in id_set:
                calls[iid] = draw_child(rand, calls[dad], calls[mom])
            else:
                # A sibship's parents aren't genotyped, but the siblings still
                # have to look like siblings, so draw one parental pair here.
                key = (variant_idx, dad)
                if key not in calls:
                    calls[key] = (sum(1 for _ in range(2) if rand.random() < freq),
                                  sum(1 for _ in range(2) if rand.random() < freq))
                pair = calls[key]
                calls[iid] = draw_child(rand, pair[0], pair[1])
        for iid in sample_ids:
            roll = rand.random()
            if roll < 0.04:
                calls[iid] = None
            elif roll < 0.07:
                calls[iid] = rand.choice([0, 1, 2])
        lines.append('1\t%d\tv%d\tA\tG\t.\t.\t.\tGT\t%s' %
                     (variant_idx + 1, variant_idx,
                      '\t'.join(geno_str(calls[iid]) for iid in sample_ids)))

    # chrX and chrMT must be dropped by --dfam.
    tail = []
    for chrom, pos, vid in (('X', 60000000, 'vx'), ('MT', 1, 'vmt')):
        calls = [rand.choice([0, 1, 2]) for _ in sample_ids]
        tail.append('%s\t%d\t%s\tA\tG\t.\t.\t.\tGT\t%s' %
                    (chrom, pos, vid,
                     '\t'.join(geno_str(c) for c in calls)))

    def write_vcf(path, flip):
        with open(path, 'w') as f:
            f.write('##fileformat=VCFv4.2\n')
            for chrom in ('1', 'X', 'MT'):
                f.write('##contig=<ID=%s>\n' % chrom)
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n')
            f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
                    '\t%s\n' % '\t'.join(sample_ids))
            for line in lines + tail:
                if not flip:
                    f.write(line + '\n')
                    continue
                g = line.split('\t')
                g[3], g[4] = g[4], g[3]
                for i in range(9, len(g)):
                    g[i] = {'0/0': '1/1', '0/1': '0/1', '1/1': '0/0',
                            './.': './.'}[g[i]]
                f.write('\t'.join(g) + '\n')

    write_vcf('dfam.vcf', False)
    write_vcf('dfam_flipped.vcf', True)


if __name__ == '__main__':
    main()

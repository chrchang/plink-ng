#!/usr/bin/env python3
"""Emit a small chrX .vcf plus a matching sex file.

Half the male calls are deliberately heterozygous, i.e. het haploid.  Those are
malformed input, and every plink command that reports per-variant counts is
expected to treat them as missing.
"""
import random
import sys

seed = int(sys.argv[1])
sample_ct = int(sys.argv[2])
variant_ct = int(sys.argv[3])
random.seed(seed)

samples = ['S%u' % (i,) for i in range(sample_ct)]
sexes = [random.choice([1, 2]) for _ in range(sample_ct)]

with open('chrx.vcf', 'w') as f:
    f.write('##fileformat=VCFv4.2\n')
    f.write('##contig=<ID=X,length=156040895>\n')
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
    f.write('\t'.join(samples))
    f.write('\n')
    for uii in range(variant_ct):
        calls = []
        for sample_idx in range(sample_ct):
            urand = random.random()
            if urand < 0.05:
                calls.append('./.')
            elif sexes[sample_idx] == 1:
                # 20% of the male calls are het haploid.
                if random.random() < 0.2:
                    calls.append('0/1')
                else:
                    calls.append(random.choice(['0/0', '1/1']))
            else:
                calls.append(random.choice(['0/0', '0/1', '1/1']))
        # Stay clear of the PAR regions so that --split-par is not needed.
        f.write('X\t%u\tx%u\tA\tG\t.\t.\t.\tGT\t' % (20000000 + uii * 1000, uii))
        f.write('\t'.join(calls))
        f.write('\n')

with open('chrx_sex.txt', 'w') as f:
    f.write('#IID\tSEX\n')
    for sample_id, sex in zip(samples, sexes):
        f.write('%s\t%u\n' % (sample_id, sex))

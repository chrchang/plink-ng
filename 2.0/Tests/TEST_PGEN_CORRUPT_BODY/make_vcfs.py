#!/usr/bin/env python3
# Writes tmp_hc.vcf (hardcalls: phased, rare, 3-allele, phased 5-allele) and
# tmp_ds.vcf (fixed-width dosages, phased dosages) for 24 samples.
N = 24
HDR = ("##fileformat=VCFv4.3\n##contig=<ID=1,length=1000>\n"
       '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
       '##FORMAT=<ID=DS,Number=A,Type=Float,Description="Dosage">\n'
       '##FORMAT=<ID=HDS,Number=.,Type=Float,Description="Haplotype dosage">\n'
       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join("s%d" % i for i in range(N)) + "\n")

def row(pos, alt, fmt, calls):
    return "1\t%d\tv%d\tA\t%s\t.\t.\t.\t%s\t%s\n" % (pos, pos, alt, fmt, "\t".join(calls))

hc = [row(100, "C", "GT", ["0|1" if i % 3 == 0 else ("1|0" if i % 3 == 1 else "0|0") for i in range(N)]),
      row(200, "C", "GT", ["0/1" if i in (2, 19) else "0/0" for i in range(N)])]
calls = []
for i in range(N):
    if i in (1, 7):
        calls.append("0/2")
    elif i == 4:
        calls.append("1/2")
    else:
        calls.append("0/1" if i % 2 else "0/0")
hc.append(row(300, "C,G", "GT", calls))
calls = []
for i in range(N):
    if i in (3, 11):
        calls.append("0|3")
    elif i == 5:
        calls.append("2|4")
    elif i == 8:
        calls.append("0/4")
    else:
        calls.append("1|0" if i % 2 else "0|0")
hc.append(row(400, "C,G,T,AA", "GT", calls))
# mostly ref/ALT1 hets with a single ref/ALT2 het: the ALT2 track is stored as
# a sample list
hc.append(row(500, "C,G", "GT", ["0/2" if i == 6 else ("0/1" if i % 4 else "0/0") for i in range(N)]))
ds = [row(500, "C", "GT:DS", ["0/1:%.3f" % (0.2 + (i % 5) * 0.3) if i % 2 else "0/0:%.3f" % (0.1 + (i % 3) * 0.1) for i in range(N)]),
      row(600, "C", "GT:DS:HDS", ["0|1:0.9:0.1,0.8" if i % 2 else "0|0:0.3:0.1,0.2" for i in range(N)])]
open("tmp_hc.vcf", "w").write(HDR + "".join(hc))
open("tmp_ds.vcf", "w").write(HDR + "".join(ds))

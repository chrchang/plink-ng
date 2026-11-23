#!/usr/bin/env python3
"""Writes a VCF for the --vcf/--bcf parity test.

usage: make_vcf.py <family> <sample_ct> <out.vcf>

Each "special" genotype is placed in its own record, at each of a set of
sample positions: the ends, and both sides of every SIMD-vector / genovec-word
boundary that fits.  All other samples in the record get unremarkable
background calls, so the special sample is the only thing in the record that
can trigger a given code path (e.g. the only phased call, or the only call
that fails a quality filter).

Families:
  geno: GT:GQ:DP, for the quality filters, --vcf-half-call, --polyploid-mode,
        mixed ploidy, multiallelic and int16 GT; also some record-level cases
        (no GT, REF N, QUAL, FILTER, >3 alleles).
  ds:   GT:DS:HDS, for dosage=DS / dosage=HDS.
  dsbad: GT:DS:HDS with a partially-missing HDS, an error under dosage=HDS.
  gp:   GT:GP, for dosage=GP and --import-dosage-certainty.
"""

import sys

ALLELES_2 = ("A", "C")
ALLELES_3 = ("A", "C,G")


def alts70():
    out = []
    for i in range(1, 70):
        a = "C"
        j = i
        while j:
            a += "ACGT"[j % 4]
            j //= 4
        out.append(a)
    return ("A", ",".join(out))


ALLELES_70 = alts70()

GENO_KINDS = [
    (ALLELES_2, "0|1:50:30"),
    (ALLELES_2, "1|0:50:30"),
    (ALLELES_2, "0|1:5:30"),     # fails --vcf-min-gq
    (ALLELES_2, "0|1:50:3"),     # fails --vcf-min-dp
    (ALLELES_2, "0|1:50:90"),    # fails --vcf-max-dp
    (ALLELES_2, "0|1:.:."),      # missing GQ/DP passes
    (ALLELES_2, "1|1:50:30"),
    (ALLELES_2, "0|0:5:30"),
    (ALLELES_2, "0/1:5:30"),
    (ALLELES_2, "1/1:50:90"),
    (ALLELES_2, "./.:50:30"),
    (ALLELES_2, ".|.:50:30"),
    (ALLELES_2, "0/.:50:30"),
    (ALLELES_2, "./1:50:30"),
    (ALLELES_2, "1|.:50:30"),
    (ALLELES_2, ".|1:5:30"),
    (ALLELES_2, "1:50:30"),      # haploid among diploids
    (ALLELES_2, "0:5:30"),
    (ALLELES_2, ".:50:30"),
    (ALLELES_2, "0/1/1:50:30"),  # triploid
    (ALLELES_2, "0|1|1:50:30"),
    (ALLELES_3, "1|2:50:30"),
    (ALLELES_3, "2|1:5:30"),
    (ALLELES_3, "2|.:50:30"),
    (ALLELES_3, "2/2:50:3"),
    (ALLELES_70, "0|69:50:30"),  # int16 GT
    (ALLELES_70, "69|1:5:30"),
    (ALLELES_70, ".|69:50:30"),
    (ALLELES_70, "69/69:50:30"),
    (ALLELES_70, "69:50:30"),
]

DS_KINDS = [
    (ALLELES_2, "0|1:0.9:0.1,0.8"),
    (ALLELES_2, "1|0:1.1:0.9,0.2"),
    (ALLELES_2, "0|1:1:0,1"),
    (ALLELES_2, "0/1:1.3:."),
    (ALLELES_2, "./.:0.7:."),
    (ALLELES_2, "./.:0.7:0.3,0.4"),
    (ALLELES_2, "0/1:.:."),
    (ALLELES_2, "1:0.8:0.8"),
    (ALLELES_2, "0:0.3:."),
    (ALLELES_2, "0/0:0.004:."),
    (ALLELES_2, "1/1:1.999:."),
    # (Multiallelic dosage import is not yet supported by either reader.)
]

# Rejected (identically) by both readers under dosage=HDS; checked as an
# expected failure.
DSBAD_KINDS = [
    (ALLELES_2, "0|1:1:0.4,."),
]

GP_KINDS = [
    (ALLELES_2, "0/1:0.1,0.8,0.1"),
    (ALLELES_2, "0/1:0.4,0.3,0.3"),
    (ALLELES_2, "0|1:0.05,0.9,0.05"),
    (ALLELES_2, "./.:0.2,0.5,0.3"),
    (ALLELES_2, "0/1:."),
    (ALLELES_2, "1:0.3,0.7"),
    (ALLELES_2, "0:0.55,0.45"),
    (ALLELES_2, "1/1:0,0.001,0.999"),
]

BG_GT = ["0/0", "0/1", "1/1", "0/0", "0/1", "0/0", "1/1"]
BG_DS = {"0/0": "0", "0/1": "1", "1/1": "2"}
BG_GP = {"0/0": "1,0,0", "0/1": "0,1,0", "1/1": "0,0,1"}


def background(family, i):
    gt = BG_GT[i % len(BG_GT)]
    if family == "geno":
        return "%s:%d:%d" % (gt, 30 + (i * 7) % 40, 15 + (i * 3) % 20)
    if family in ("ds", "dsbad"):
        return "%s:%s:." % (gt, BG_DS[gt])
    return "%s:%s" % (gt, BG_GP[gt])


def positions(n):
    pos = {0, 1, 2, n // 2, n - 2, n - 1}
    for k in (4, 8, 16, 32, 64, 128):
        pos.update((k - 1, k, k + 1))
    return sorted(p for p in pos if 0 <= p < n)


HEADERS = {
    "geno": ['##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">',
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
             '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">'],
    "ds": ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
           '##FORMAT=<ID=DS,Number=A,Type=Float,Description="Dosage">',
           '##FORMAT=<ID=HDS,Number=.,Type=Float,Description="Haplotype dosage">'],
    "gp": ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
           '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype probabilities">'],
}
FORMATS = {"geno": "GT:GQ:DP", "ds": "GT:DS:HDS", "gp": "GT:GP"}
KINDS = {"geno": GENO_KINDS, "ds": DS_KINDS, "dsbad": DSBAD_KINDS, "gp": GP_KINDS}
HEADERS["dsbad"] = HEADERS["ds"]
FORMATS["dsbad"] = FORMATS["ds"]


def main():
    family, n, out_fname = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    bg = [background(family, i) for i in range(n)]
    fmt = FORMATS[family]
    rows = []

    def add(alleles, samples, fmt=fmt, qual=".", filt=".", ref=None):
        pos = 100 * (len(rows) + 1)
        rows.append("\t".join(["1", str(pos), "v%d" % (len(rows) + 1),
                               ref or alleles[0], alleles[1], qual, filt, ".",
                               fmt] + samples))

    add(ALLELES_2, bg)
    for kind_idx, (alleles, special) in enumerate(KINDS[family]):
        for p in positions(n):
            samples = list(bg)
            samples[p] = special
            add(alleles, samples)
    if family == "geno":
        # Record-level cases, with a phased het in the first sample so the
        # GT scan has something to find.
        ph0 = ["0|1:50:30"] + bg[1:]
        add(ALLELES_2, ["50:30"] * n, fmt="GQ:DP")       # no GT
        add(ALLELES_2, ["50"] * n, fmt="GQ")
        add(ALLELES_2, ph0, ref="N")                     # REF N
        add(ALLELES_2, ph0, qual="5")
        add(ALLELES_2, ph0, qual="50")
        add(ALLELES_2, ph0, filt="PASS")
        add(("A", "C,G,T"), ["1|3:50:30"] + bg[1:])      # 4 alleles
    with open(out_fname, "w") as f:
        f.write("##fileformat=VCFv4.3\n##contig=<ID=1,length=100000000>\n")
        f.write("\n".join(HEADERS[family]) + "\n")
        f.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                           "FILTER", "INFO", "FORMAT"]
                          + ["s%d" % (i + 1) for i in range(n)]) + "\n")
        f.write("\n".join(rows) + "\n")


if __name__ == "__main__":
    main()

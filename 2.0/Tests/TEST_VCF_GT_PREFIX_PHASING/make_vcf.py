#!/usr/bin/env python3
"""Writes VCFv4.4 GT-prefix-phasing fixtures for run_tests.sh.

usage: make_vcf.py <family> <sample count> <prefixed out.vcf> <canonical out.vcf>

VCFv4.4+ lets each GT value start with an optional phasing indicator, e.g.
"|0/1" or "/0|1"; plink2 skips it, so for a diploid call phase comes from
the indicator between the two alleles.  The canonical file is the prefixed
one with every leading indicator removed, so the two must import
identically.

Each record cycles through a list of genotypes with a per-record offset, so
every genotype lands at the first and last sample and at a variety of
positions relative to genovec-word boundaries.  Records whose only phased het
is prefixed check the phased-het pre-scans; the >10-allele records cover the
long-allele parse and pre-scan paths.  Families:
  geno: GT:GQ:DP (GQ/DP give the quality-filter tests something to filter)
  gt:   GT only, so the 32-sample SWAR fast path is attempted
  ds:   GT:DS, where a prefix must not make a diploid call look haploid
        (haploid DS values are doubled)
  hds:  GT:DS:HDS
  gp:   GT:GP, including diploid calls like 0.1,0.1,0.8 that a prefix-induced
        haploid misreading would wrongly force-missing under
        --import-dosage-certainty
  poly: a real triploid call, which must still fail without
        --polyploid-mode missing
"""

import sys

ALTS12 = ["C", "G", "T", "AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC"]

DIPLOID = ["|0/0", "/0|1", "|0|1", "|1/0", "/1|0", "|0/1", "/1/1", "|./.",
           "/.|.", "0/1", "0|1", "1|1", "0/0"]
HALF = ["|./1", "/.|1", "|1/.", ".|1"]
HAPLOID = ["|1", "/0", "|.", "1", "0"]
MULTI3 = ["|0/2", "/2|1", "|1|2", "/2/2", "0|2", "|2|0"]
MULTI12 = ["|11/3", "/2|11", "|10|0", "/11/11", "|0/10", "11|1", "|.|11"]


def strip(gt):
    return gt[1:] if gt[0] in "/|" else gt


def cyc(vals, n, offset):
    return [vals[(i + offset) % len(vals)] for i in range(n)]


def main():
    fam, n, out_pref, out_canon = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
    # (alts, genotype list or explicit per-sample list)
    recs = []
    for off in range(3):
        recs.append((["C"], cyc(DIPLOID, n, off)))
    recs.append((["C"], cyc(DIPLOID + HALF, n, 5)))
    if fam in ("geno", "gt"):
        # first sample is a prefixed half-call, so the default
        # --vcf-half-call error mode always fails
        recs.append((["C"], ["|./1"] + cyc(DIPLOID, n - 1, 4)))
    recs.append((["C"], cyc(DIPLOID + HAPLOID, n, 1)))
    recs.append((["C"], cyc(HAPLOID, n, 2)))
    # only phased het is prefixed: "/0|1" in one, "|0|1" in another, at the
    # last sample.  And one where "|0/1" is the only '|' (must stay unphased).
    for special in ("/0|1", "|0|1", "|1/0"):
        recs.append((["C"], ["0/0"] * (n - 1) + [special]))
    if fam not in ("ds", "hds", "gp"):
        recs.append((["C", "G"], cyc(MULTI3 + DIPLOID, n, 0)))
        recs.append((["C", "G"], ["0/0"] * (n - 1) + ["/2|1"]))
        recs.append((ALTS12, cyc(MULTI12 + DIPLOID, n, 0)))
        recs.append((ALTS12, ["0/0"] * (n - 1) + ["/2|11"]))
        recs.append((ALTS12, ["0/0"] * (n - 1) + ["|11|11"]))
        recs.append((ALTS12, cyc(MULTI12 + HAPLOID, n, 3)))
    if fam == "poly":
        recs.append((["C"], ["0/0"] * (n - 1) + ["|0/1/1"]))

    fmt = {"geno": "GT:GQ:DP", "ds": "GT:DS", "hds": "GT:DS:HDS", "gp": "GT:GP"}.get(fam, "GT")
    gp_dip = [(0.1, 0.1, 0.8), (0.8, 0.1, 0.1), (0.3, 0.4, 0.3), (0.05, 0.9, 0.05), (0.2, 0.2, 0.6)]
    gp_hap = [(0.1, 0.9), (0.6, 0.4), (0.5, 0.5), (0.95, 0.05)]

    def extra(i, j, haploid):
        a = ((i * 7 + j * 3) % 9) / 10.0 + 0.05
        b = ((i * 5 + j * 11) % 9) / 10.0 + 0.05
        if fam == "geno":
            return ":%d:%d" % ((i * 7 + j * 3) % 60, (i * 5 + j * 11) % 50)
        if fam == "ds":
            return ":%g" % (a if haploid else a + b)
        if fam == "hds":
            return (":%g:%g" % (a, a)) if haploid else (":%g:%g,%g" % (a + b, a, b))
        if fam == "gp":
            vals = gp_hap[(i + j) % len(gp_hap)] if haploid else gp_dip[(i * 3 + j) % len(gp_dip)]
            return ":" + ",".join("%g" % x for x in vals)
        return ""

    header = [
        "##fileformat=VCFv4.4",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">',
        '##FORMAT=<ID=DS,Number=A,Type=Float,Description="DS">',
        '##FORMAT=<ID=HDS,Number=.,Type=Float,Description="HDS">',
        '##FORMAT=<ID=GP,Number=G,Type=Float,Description="GP">',
        "##contig=<ID=1,length=1000000>",
        "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
                  + ["s%d" % (j + 1) for j in range(n)]),
    ]
    with open(out_pref, "w") as fp, open(out_canon, "w") as fc:
        fp.write("\n".join(header) + "\n")
        fc.write("\n".join(header) + "\n")
        for i, (alts, gts) in enumerate(recs):
            gts = gts[:n]
            def field(gt, j, canon):
                haploid = ("/" not in strip(gt)) and ("|" not in strip(gt))
                return (strip(gt) if canon else gt) + extra(i, j, haploid)

            cols = ["1", str(100 * (i + 1)), "v%d" % (i + 1), "A", ",".join(alts), ".", ".", ".", fmt]
            fp.write("\t".join(cols + [field(g, j, False) for j, g in enumerate(gts)]) + "\n")
            fc.write("\t".join(cols + [field(g, j, True) for j, g in enumerate(gts)]) + "\n")


main()

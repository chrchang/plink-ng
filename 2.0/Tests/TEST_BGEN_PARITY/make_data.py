#!/usr/bin/env python3
"""Writes matched VCF and BGEN fixtures for the --bgen/--vcf parity test.

usage: make_data.py <family> <sample_ct> <bits|L1> <out_prefix> <encodings>

  family:    gp  (unphased; the VCF carries GP, imported with dosage=GP)
             hds (phased; the VCF carries HDS, imported with dosage=HDS)
  bits:      BGEN layout-2 bit depth 1..28, or L1 for layout 1 (gp only)
  encodings: comma-separated subset of
               zlib, none, zstd  (layout-2 compression, embedded sample IDs)
               noids             (zlib, no sample ID block; writes .sample)
               reflast           (zlib, alleles in reverse order and
                                  probabilities mirrored, for 'ref-last')
               l1zlib, l1none    (layout 1)

Every probability is first quantized to the BGEN scale (k / (2^bits - 1), or
k / 32768 for layout 1), and the VCF gets the decimal value of exactly that
rational, so both files encode the same numbers.  Each special genotype is
placed alone at the first and last sample and on both sides of every
SIMD-vector / genovec-word boundary that fits; all other samples in the
record get hardcall-valued background probabilities (which keeps the 8- and
16-bit all-hardcall fast paths live until the special sample).  Specials
that change ploidy are placed both in a diploid background (min ploidy !=
max ploidy) and in a background of their own ploidy (fixed-ploidy path).
"""

import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
import bgen_writer  # noqa: E402

F = Fraction

# Each special: (ploidy, missing, probabilities).  Unphased diploid
# probabilities are (P(0/0), P(0/1), P(1/1)); unphased haploid (P(0), P(1));
# phased: per-haplotype ALT dosages.  Ploidy 0 / missing carry no values.
GP_SPECIALS = [
    (2, False, (F(2, 10), F(7, 10), F(1, 10))),
    (2, False, (F(0), F(95, 100), F(5, 100))),
    (2, False, (F(45, 100), F(1, 10), F(45, 100))),
    (2, False, (F(1, 3), F(1, 3), F(1, 3))),
    (2, False, (F(9, 10), F(1, 10), F(0))),
    (2, False, (F(6, 10), F(2, 10), F(2, 10))),   # --import-dosage-certainty 0.6 tie at 4 bits
    (2, False, (F(5, 100), F(5, 100), F(9, 10))),
    (2, False, (F(0), F(1), F(0))),
    (2, False, (F(1, 15), F(13, 15), F(1, 15))),
    (2, True, None),
    (1, False, (F(3, 10), F(7, 10))),
    (1, False, (F(6, 10), F(4, 10))),
    (1, False, (F(0), F(1))),
    (1, True, None),
    (0, False, None),
]

# Layout 1 stores all three probabilities, which need not sum to 1; plink2
# normalizes by their sum.  (Layout 1 only.)
GP_L1_UNNORMALIZED = [
    (2, False, (F(3, 10), F(5, 10), F(1, 10))),
    (2, False, (F(2, 10), F(2, 10), F(2, 10))),
    (2, False, (F(1, 100), F(1, 100), F(1, 100))),
    (2, False, (F(0), F(3, 10), F(6, 10))),
]

HDS_SPECIALS = [
    (2, False, (F(0), F(1))),
    (2, False, (F(1), F(0))),
    (2, False, (F(9, 10), F(2, 10))),
    (2, False, (F(1, 2), F(1, 2))),
    (2, False, (F(95, 100), F(5, 100))),
    (2, False, (F(2, 100), F(3, 100))),
    (2, False, (F(6, 10), F(6, 10))),
    (2, False, (F(1, 15), F(14, 15))),
    (2, True, None),
    (1, False, (F(8, 10),)),
    (1, False, (F(4, 10),)),
    (1, False, (F(1, 2),)),       # haploid near 0.5 (--dosage-erase-threshold)
    (1, False, (F(1),)),
    (1, True, None),
    (0, False, None),
]

GP_BG = [(F(1), F(0), F(0)), (F(0), F(1), F(0)), (F(0), F(0), F(1)), (F(1), F(0), F(0))]
GP_BG_HAP = [(F(1), F(0)), (F(0), F(1))]
HDS_BG = [(F(0), F(0)), (F(0), F(1)), (F(1), F(0)), (F(1), F(1)), (F(0), F(0))]
HDS_BG_HAP = [(F(0),), (F(1),)]


def positions(n):
    pos = {0, 1, 2, n // 2, n - 2, n - 1}
    for k in (4, 8, 16, 32, 64):
        pos.update((k - 1, k, k + 1))
    return sorted(p for p in pos if 0 <= p < n)


def background(family, ploidy, i):
    if family == "gp":
        bg = GP_BG if ploidy == 2 else GP_BG_HAP
    else:
        bg = HDS_BG if ploidy == 2 else HDS_BG_HAP
    return (ploidy, False, bg[i % len(bg)])


def records(family, n, diploid_only):
    specials = GP_SPECIALS if family == "gp" else HDS_SPECIALS
    if diploid_only:
        specials = [sp for sp in specials if sp[0] == 2] + GP_L1_UNNORMALIZED
    out = [[background(family, 2, i) for i in range(n)]]
    for sp in specials:
        bg_ploidies = [2] if sp[0] == 2 else [2, 1]
        for bgp in bg_ploidies:
            for p in positions(n):
                rec = [background(family, bgp, i) for i in range(n)]
                rec[p] = sp
                out.append(rec)
    return out


def quantize(probs, denom):
    """Round cumulative sums, so the numerators are nonnegative and sum to
    denom."""
    ns = []
    prev = 0
    cum = F(0)
    for p in probs[:-1]:
        cum += p
        c = int(cum * denom + F(1, 2))
        ns.append(c - prev)
        prev = c
    ns.append(denom - prev)
    return ns


def dec(num, denom):
    if num == 0:
        return "0"
    if num == denom:
        return "1"
    return repr(num / denom)


def build(family, n, bits):
    """Returns (vcf_lines, variants_layout2_or_1, layout)."""
    layout1 = bits == "L1"
    denom = 32768 if layout1 else (1 << bits) - 1
    recs = records(family, n, layout1)
    vcf_rows = []
    variants = []
    for ridx, rec in enumerate(recs):
        vid = "v%d" % (ridx + 1)
        pos = 100 * (ridx + 1)
        cells = []
        numers = []
        l1probs = []
        ploidies = []
        missings = []
        for ploidy, missing, probs in rec:
            ploidies.append(ploidy)
            missings.append(missing)
            if family == "gp":
                if ploidy == 0 or missing:
                    cells.append(("./." if ploidy == 2 else ".") + ":.")
                    numers.append([0] * ploidy)
                    l1probs.append((0, 0, 0))
                    continue
                if layout1 and sum(probs) != 1:
                    # Round independently; keep the sum odd so the
                    # normalized dosage can't land exactly on a rounding tie.
                    q = [int(p * denom + F(1, 2)) for p in probs]
                    if sum(q) % 2 == 0:
                        q[0] += 1
                else:
                    q = quantize(probs, denom)
                    if layout1 and (q[1] % 2):
                        # With a 32768 sum, odd P(0/1) numerators put the
                        # dosage exactly halfway between two 1/16384 steps.
                        # Layout-1 import rounds those ties to even, while
                        # VCF GP import rounds them up (by design); avoid
                        # them.
                        q[1] -= 1
                        q[0] += 1
                gt = "0/0" if ploidy == 2 else "0"
                cells.append(gt + ":" + ",".join(dec(x, denom) for x in q))
                numers.append(q[:-1])
                if layout1:
                    l1probs.append(tuple(q))
            else:
                if ploidy == 0 or missing:
                    cells.append(("./." if ploidy == 2 else ".") + ":.")
                    numers.append([0] * ploidy)
                    continue
                # stored value is P(haplotype carries REF)
                ref_q = [int((1 - h) * denom + F(1, 2)) for h in probs]
                alt = [denom - r for r in ref_q]
                gt = "|".join("1" if 2 * a > denom else "0" for a in alt)
                cells.append(gt + ":" + ",".join(dec(a, denom) for a in alt))
                numers.append(ref_q)
        vcf_rows.append("\t".join(["1", str(pos), vid, "A", "C", ".", ".", ".",
                                   "GT:GP" if family == "gp" else "GT:HDS"] + cells))
        v = {"chrom": "1", "pos": pos, "varid": vid, "rsid": vid,
             "alleles": ["A", "C"]}
        if layout1:
            assert all(p == 2 for p in ploidies)
            v["probs"] = l1probs
        else:
            v.update(phased=(family == "hds"), ploidy=ploidies,
                     missing=missings, numers=numers)
        variants.append(v)
    return vcf_rows, variants, denom


def mirror(v, family, denom):
    """Same genotypes with the allele order reversed (for 'ref-last')."""
    w = dict(v)
    w["alleles"] = list(reversed(v["alleles"]))
    new_numers = []
    for ploidy, missing, ns in zip(v["ploidy"], v["missing"], v["numers"]):
        if ploidy == 0 or missing:
            new_numers.append(list(ns))
        elif family == "hds":
            new_numers.append([denom - x for x in ns])
        elif ploidy == 1:
            new_numers.append([denom - ns[0]])
        else:
            n2 = denom - ns[0] - ns[1]
            new_numers.append([n2, ns[1]])
    w["numers"] = new_numers
    return w


def main():
    family, n, bits, prefix, encodings = sys.argv[1:6]
    n = int(n)
    if bits != "L1":
        bits = int(bits)
    encodings = encodings.split(",")
    sample_ids = ["s%d" % (i + 1) for i in range(n)]
    vcf_rows, variants, denom = build(family, n, bits)
    with open(prefix + ".vcf", "w") as f:
        f.write("##fileformat=VCFv4.3\n##contig=<ID=1,length=100000000>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        if family == "gp":
            f.write('##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype probabilities">\n')
        else:
            f.write('##FORMAT=<ID=HDS,Number=.,Type=Float,Description="Haplotype dosages">\n')
        f.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                           "INFO", "FORMAT"] + sample_ids) + "\n")
        f.write("\n".join(vcf_rows) + "\n")
    for enc in encodings:
        out = "%s.%s.bgen" % (prefix, enc)
        if enc in ("zlib", "none", "zstd"):
            if enc == "zstd" and not bgen_writer.zstd_available():
                continue
            bgen_writer.write_bgen(out, sample_ids, variants, 2, bits, enc)
        elif enc == "noids":
            bgen_writer.write_bgen(out, sample_ids, variants, 2, bits, "zlib",
                                   embed_ids=False)
            bgen_writer.write_sample_file(prefix + ".noids.sample", sample_ids)
        elif enc == "reflast":
            bgen_writer.write_bgen(out, sample_ids,
                                   [mirror(v, family, denom) for v in variants],
                                   2, bits, "zlib")
        elif enc in ("l1zlib", "l1none"):
            bgen_writer.write_bgen(out, sample_ids, variants, 1, 16,
                                   "zlib" if enc == "l1zlib" else "none")
        else:
            raise ValueError(enc)


if __name__ == "__main__":
    main()

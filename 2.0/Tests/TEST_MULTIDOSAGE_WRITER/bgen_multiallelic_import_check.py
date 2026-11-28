#!/usr/bin/env python3
"""Independent check of plink2 --bgen import of variants with more than two
alleles.

Usage:
  bgen_multiallelic_import_check.py <.bgen> <.pvar> <dumpcalls output>
      <ref-first | ref-last> <hct> <erase> <certainty> [max_allele_ct]

<hct> and <erase> are the --hard-call-threshold and --dosage-erase-threshold
distances in 16384ths of a copy (plink2's integer values; the import default
hct is 1638, erase 0).  <certainty> is the --import-dosage-certainty value as
given on the command line (0 if unused).  The .bgen variant IDs must be the plink2 variant IDs (rsid
field).  Each variant's .pvar alleles must be the .bgen alleles in the
expected order (REF moved from the end with ref-last).  Biallelic variants'
genotypes are skipped.  With [max_allele_ct]
(--import-max-alleles), variants with more alleles must be absent.

Expected stored record for each sample, per the plink2 definitions:
  * missing flag, or ploidy 0: missing hardcall, no dosage.
  * --import-dosage-certainty: missing unless the largest genotype probability
    (haploid: allele probability) is at least the threshold; for phased
    diploid samples, the genotype distribution with maximal heterozygosity
    (one haplotype takes the highest-dosage alleles up to one copy, the other
    the rest, independently) is used.
  * allele dosages u: unphased diploid, sum over genotypes of P(g) * (copies
    of the allele in g), scaled to 32768 per sample; haploid, the allele
    probabilities scaled to 16384, then doubled; phased diploid, each
    haplotype scaled to 16384, u = h1 + h2 and d = h1 - h2.  Each scaling uses
    largest-remainder rounding, ties to the lower output allele index.
  * hardcall: each allele's dosage rounded to the nearest copy (ties up); it
    is missing unless the counts sum to 2 and half the summed absolute
    rounding error is at most hct.  A heterozygous call a < b of a phased
    sample with a nonzero haplotype split is phased iff
    |(d_b - d_a)| > 16384, with b on the first haplotype iff positive.
  * the dosage (and dosage-phase, stored iff the split is nonzero) is erased
    when it is within the erase threshold: without a split, the summed
    absolute rounding error is at most 2 * erase; with one, each haplotype's
    summed absolute rounding error (32768 = one copy) is at most 2 * erase.
"""

import shutil
import struct
import subprocess
import sys
import zlib

try:
    from compression import zstd as _zstd  # Python 3.14+
except ImportError:
    _zstd = None


def fail(msg):
    sys.exit("bgen_multiallelic_import_check: " + msg)


def zstd_decompress(data):
    if _zstd is not None:
        return _zstd.decompress(data)
    if shutil.which("zstd") is None:
        return None
    return subprocess.run(["zstd", "-q", "-d", "-c"], input=data,
                          stdout=subprocess.PIPE, check=True).stdout


def unpack_bits(data, bits, ct):
    mask = (1 << bits) - 1
    out = []
    acc = 0
    acc_bits = 0
    byte_iter = iter(data)
    for _ in range(ct):
        while acc_bits < bits:
            acc |= next(byte_iter) << acc_bits
            acc_bits += 8
        out.append(acc & mask)
        acc >>= bits
        acc_bits -= bits
    return out


def lr_round(numers, total, target):
    # largest remainder; ties to the lower index
    q = [(n * target) // total for n in numers]
    rems = [(n * target) % total for n in numers]
    deficit = target - sum(q)
    order = sorted(range(len(numers)), key=lambda i: (-rems[i], i))
    for i in order[:deficit]:
        q[i] += 1
    return q


def greedy_max_prob(sums, max_val):
    # sums: per-allele haplotype-numerator sums (2 * max_val total)
    xs = sorted([x for x in sums if x], reverse=True)
    h1 = []
    h2 = []
    cap = max_val
    for x in xs:
        take = min(x, cap)
        h1.append(take)
        h2.append(x - take)
        cap -= take
    best = 0
    for i in range(len(xs)):
        best = max(best, h1[i] * h2[i])
        for j in range(i + 1, len(xs)):
            best = max(best, h1[i] * h2[j] + h1[j] * h2[i])
    return best  # times max_val^2


def expected_token(u, d, hct, erase):
    k = len(u)
    r = [(x + 8192) // 16384 for x in u]
    err = sum(abs(u[i] - 16384 * r[i]) for i in range(k))
    if (sum(r) != 2) or (err > 2 * hct):
        call = "."
    else:
        alleles = [i for i in range(k) for _ in range(r[i])]
        a, b = alleles
        call = "%d/%d" % (a, b)
        if (d is not None) and (a != b):
            diff = d[b] - d[a]
            if abs(diff) > 16384:
                call = ("%d|%d" % (b, a)) if diff > 0 else ("%d|%d" % (a, b))
    if d is None:
        erased = err <= 2 * erase
    else:
        erased = True
        for sign in (1, -1):
            hap = [u[i] + sign * d[i] for i in range(k)]
            dist = sum(abs(h - 32768 * ((h + 16384) // 32768)) for h in hap)
            erased = erased and (dist <= 2 * erase)
    if erased:
        return call
    tok = call + "=" + ",".join(str(x) for x in u)
    if d is not None:
        tok += "~" + ",".join(str(x) for x in d)
    return tok


def main():
    bgen_fname, pvar_fname, calls_fname, ref_mode, hct, erase, certainty = sys.argv[1:8]
    pvar_alleles = {}
    for line in open(pvar_fname):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        pvar_alleles[f[2]] = [f[3]] + f[4].split(",")
    hct = int(hct)
    erase = int(erase)
    # plink2 shrinks the --import-dosage-certainty value by 2^-44 relative
    certainty = float(certainty) * (1.0 - 2.0 ** -44)
    max_allele_ct = int(sys.argv[8]) if len(sys.argv) > 8 else 65535
    filtered_ct = 0
    ref_last = (ref_mode == "ref-last")
    calls = {}
    for line in open(calls_fname):
        f = line.rstrip("\n").split("\t")
        calls[f[0]] = f[2:]
    raw = open(bgen_fname, "rb").read()
    offset, header_len, variant_ct, sample_ct = struct.unpack_from("<IIII", raw, 0)
    flags = struct.unpack_from("<I", raw, 4 + header_len - 4)[0]
    compression = flags & 3
    if ((flags >> 2) & 15) != 2:
        fail("expected layout 2")
    pos = offset + 4
    stats = {"variants": 0, "samples": 0, "dosages": 0, "dphases": 0, "phased": 0, "missing": 0, "uncertain": 0, "haploid": 0}
    for _ in range(variant_ct):
        (slen,) = struct.unpack_from("<H", raw, pos)
        pos += 2 + slen
        (slen,) = struct.unpack_from("<H", raw, pos)
        varid = raw[pos + 2:pos + 2 + slen].decode()
        pos += 2 + slen
        (slen,) = struct.unpack_from("<H", raw, pos)
        pos += 2 + slen + 4
        (allele_ct,) = struct.unpack_from("<H", raw, pos)
        pos += 2
        bgen_alleles = []
        for _ in range(allele_ct):
            (alen,) = struct.unpack_from("<I", raw, pos)
            bgen_alleles.append(raw[pos + 4:pos + 4 + alen].decode())
            pos += 4 + alen
        (block_len,) = struct.unpack_from("<I", raw, pos)
        if compression == 1:
            block = zlib.decompress(raw[pos + 8:pos + 4 + block_len])
        elif compression == 0:
            block = raw[pos + 4:pos + 4 + block_len]
        else:
            block = zstd_decompress(raw[pos + 8:pos + 4 + block_len])
            if block is None:
                print("bgen_multiallelic_import_check: zstd unavailable, skipping " + bgen_fname)
                return
        pos += 4 + block_len
        if allele_ct > max_allele_ct:
            if varid in calls:
                fail(varid + " should have been filtered out")
            filtered_ct += 1
            continue
        expected_alleles = (bgen_alleles[-1:] + bgen_alleles[:-1]) if ref_last else bgen_alleles
        if pvar_alleles.get(varid) != expected_alleles:
            fail(varid + ": .pvar alleles %s, expected %s" % (pvar_alleles.get(varid), expected_alleles))
        if allele_ct == 2:
            continue
        if varid not in calls:
            fail(varid + " missing from dumpcalls output")
        got = calls[varid]
        n, k = struct.unpack_from("<IH", block, 0)
        if (n != sample_ct) or (k != allele_ct) or (len(got) != sample_ct):
            fail(varid + ": count mismatch")
        ploidy_bytes = block[8:8 + sample_ct]
        phased, bits = block[8 + sample_ct], block[9 + sample_ct]
        max_val = (1 << bits) - 1
        thresh = (1 + int(certainty * max_val)) if certainty > 0 else 0
        per_ploidy = {0: 0, 1: k - 1, 2: (2 * (k - 1)) if phased else (k * (k + 1) // 2 - 1)}
        tot = sum(per_ploidy[p & 63] for p in ploidy_bytes)
        vals = unpack_bits(block[10 + sample_ct:], bits, tot)
        out_idx = [((i + 1) if i != k - 1 else 0) if ref_last else i for i in range(k)]
        vpos = 0
        for sidx in range(sample_ct):
            p = ploidy_bytes[sidx]
            ploidy = p & 63
            cur = vals[vpos:vpos + per_ploidy[ploidy]]
            vpos += per_ploidy[ploidy]
            if (p & 128) or (not ploidy):
                expected = "."
                stats["missing"] += 1
            else:
                d = None
                if phased or (ploidy == 1):
                    haps = []
                    for h in range(ploidy):
                        hv = cur[h * (k - 1):(h + 1) * (k - 1)]
                        hv = hv + [max_val - sum(hv)]
                        if hv[-1] < 0:
                            fail(varid + ": numerator overflow in input")
                        out = [0] * k
                        for i in range(k):
                            out[out_idx[i]] = hv[i]
                        haps.append(out)
                    if ploidy == 1:
                        certain = max(haps[0]) >= thresh
                        u = [2 * x for x in lr_round(haps[0], max_val, 16384)]
                        stats["haploid"] += 1
                    else:
                        sums = [haps[0][i] + haps[1][i] for i in range(k)]
                        certain = greedy_max_prob(sums, max_val) >= thresh * max_val
                        h1 = lr_round(haps[0], max_val, 16384)
                        h2 = lr_round(haps[1], max_val, 16384)
                        u = [h1[i] + h2[i] for i in range(k)]
                        if h1 != h2:
                            d = [h1[i] - h2[i] for i in range(k)]
                else:
                    gp = cur + [max_val - sum(cur)]
                    if gp[-1] < 0:
                        fail(varid + ": numerator overflow in input")
                    certain = max(gp) >= thresh
                    numers = [0] * k
                    gidx = 0
                    for b in range(k):
                        for a in range(b + 1):
                            numers[out_idx[a]] += gp[gidx]
                            numers[out_idx[b]] += gp[gidx]
                            gidx += 1
                    u = lr_round(numers, 2 * max_val, 32768)
                if not certain:
                    expected = "."
                    stats["uncertain"] += 1
                else:
                    expected = expected_token(u, d, hct, erase)
            if got[sidx] != expected:
                fail("%s sample %d: got %s, expected %s" % (varid, sidx, got[sidx], expected))
            stats["samples"] += 1
            stats["dosages"] += ("=" in expected)
            stats["dphases"] += ("~" in expected)
            stats["phased"] += ("|" in expected)
        stats["variants"] += 1
    if (max_allele_ct < 65535) and (not filtered_ct):
        fail("no variants filtered by max_allele_ct")
    if not stats["variants"]:
        fail("no multiallelic variants checked")
    print("bgen_multiallelic_import_check: %s ok (%s)" % (bgen_fname, ", ".join("%d %s" % (stats[key], key) for key in ("variants", "samples", "dosages", "dphases", "phased", "missing", "uncertain", "haploid"))))


main()

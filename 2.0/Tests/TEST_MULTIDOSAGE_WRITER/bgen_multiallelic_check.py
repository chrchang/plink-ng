#!/usr/bin/env python3
"""Independent check of plink2 --export bgen-1.2/1.3 output for variants with
more than two alleles.

Usage:
  bgen_multiallelic_check.py <.bgen> <.pvar> <.psam> <dumpvals output>
                             <ref-first | ref-last>

The dumpvals file (multidosage_writer_test dumpvals) gives each sample's
per-allele dosages u and, when the haplotype split is known, per-allele
haplotype differences d, in .pvar allele order.  Expected bgen layout-2
contents, per the plink2_export.cc definition:
  * phased format (one probability vector per haplotype);
  * ploidy 2 on diploid chromosomes; on chrX, males get ploidy 1 unless the
    split is known; on chrY, missing females get ploidy 0, others ploidy 1
    unless the split is known; chrM likewise (ploidy 1 unless split);
  * haplotype values u + d / u - d when split, u for both haplotypes (or the
    single one) otherwise, out of 32768;
  * each haplotype's values scaled to 2^bits - 1 with largest-remainder
    rounding (ties to the lower output allele index), last allele omitted;
  * missing samples flagged, with zero probabilities.
Output allele order is taken from the .bgen variant header and matched to
the .pvar alleles by string; it must be the .pvar order (ref-first), or that
order with REF moved to the end (ref-last), as specified.  Biallelic
variants' genotype blocks are skipped.
The 'missing' column of the .sample file next to the .bgen must equal each
sample's fraction of missing calls over all variants (chrY calls of females
excluded).
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
    sys.exit("bgen_multiallelic_check: " + msg)


def zstd_decompress(data):
    if _zstd is not None:
        return _zstd.decompress(data)
    if shutil.which("zstd") is None:
        return None
    return subprocess.run(["zstd", "-q", "-d", "-c"], input=data,
                          stdout=subprocess.PIPE, check=True).stdout


def read_pvar(fname):
    chroms = {}
    alleles = {}
    for line in open(fname):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        chroms[f[2]] = f[0]
        alleles[f[2]] = [f[3]] + f[4].split(",")
    return chroms, alleles


def read_psam(fname):
    sexes = []
    sex_col = None
    for line in open(fname):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#"):
            sex_col = f.index("SEX") if "SEX" in f else None
            continue
        sexes.append(f[sex_col] if sex_col is not None else "NA")
    return sexes


def read_vals(fname):
    vals = {}
    for line in open(fname):
        f = line.rstrip("\n").split("\t")
        vals[f[0]] = f[2:]
    return vals


def unpack_bits(data, bits, ct):
    # little-endian bit order, streamed (whole-block big integers are
    # quadratic here)
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


def round_hap(hap, max_val):
    # hap: values summing to 32768
    numers = [(v * max_val) // 32768 for v in hap]
    rems = [(v * max_val) % 32768 for v in hap]
    deficit = max_val - sum(numers)
    order = sorted(range(len(hap)), key=lambda i: (-rems[i], i))
    for i in order[:deficit]:
        numers[i] += 1
    return numers


def main():
    bgen_fname, pvar_fname, psam_fname, vals_fname, ref_mode = sys.argv[1:6]
    if ref_mode not in ("ref-first", "ref-last"):
        fail("invalid REF mode")
    chroms, pvar_alleles = read_pvar(pvar_fname)
    sexes = read_psam(psam_fname)
    vals = read_vals(vals_fname)
    raw = open(bgen_fname, "rb").read()
    offset, header_len, variant_ct, sample_ct = struct.unpack_from("<IIII", raw, 0)
    flags = struct.unpack_from("<I", raw, header_len)[0]
    compression = flags & 3
    layout = (flags >> 2) & 15
    if layout != 2:
        fail("expected layout 2, got %d" % layout)
    if sample_ct != len(sexes):
        fail("sample count mismatch")
    any_male = "1" in sexes
    any_female = "2" in sexes
    pos = offset + 4
    checked_variants = 0
    checked_samples = 0
    ref_last_seen = set()
    missing_cts = [0] * sample_ct
    split_samples = 0
    haploid_samples = 0
    for _ in range(variant_ct):
        # plink2 leaves the first ID empty and stores the variant ID as the
        # rsid
        (slen,) = struct.unpack_from("<H", raw, pos)
        pos += 2 + slen
        (slen,) = struct.unpack_from("<H", raw, pos)
        varid = raw[pos + 2:pos + 2 + slen].decode()
        pos += 2 + slen
        (slen,) = struct.unpack_from("<H", raw, pos)
        pos += 2 + slen  # chrom
        pos += 4  # bp
        (allele_ct,) = struct.unpack_from("<H", raw, pos)
        pos += 2
        out_alleles = []
        for _ in range(allele_ct):
            (alen,) = struct.unpack_from("<I", raw, pos)
            out_alleles.append(raw[pos + 4:pos + 4 + alen].decode())
            pos += 4 + alen
        (block_len,) = struct.unpack_from("<I", raw, pos)
        if compression:
            (uncompressed_len,) = struct.unpack_from("<I", raw, pos + 4)
            comp = raw[pos + 8:pos + 4 + block_len]
            if compression == 1:
                block = zlib.decompress(comp)
            else:
                block = zstd_decompress(comp)
                if block is None:
                    print("bgen_multiallelic_check: zstd unavailable, skipping " + bgen_fname)
                    return
            if len(block) != uncompressed_len:
                fail(varid + ": uncompressed length mismatch")
        else:
            block = raw[pos + 4:pos + 4 + block_len]
        pos += 4 + block_len
        in_alleles = pvar_alleles[varid]
        if out_alleles == in_alleles:
            ref_last_seen.add(0)
        elif out_alleles == in_alleles[1:] + in_alleles[:1]:
            ref_last_seen.add(1)
        else:
            fail(varid + ": unexpected output allele order")
        y_variant = chroms[varid] in ("Y", "24", "chrY")
        for sidx, tok in enumerate(vals[varid]):
            if (tok == "M") and not (y_variant and sexes[sidx] == "2"):
                missing_cts[sidx] += 1
        if allele_ct == 2:
            continue
        out_to_in = [in_alleles.index(a) for a in out_alleles]
        n, k, pmin, pmax = struct.unpack_from("<IHBB", block, 0)
        if (n != sample_ct) or (k != allele_ct):
            fail(varid + ": block header mismatch")
        ploidy_bytes = block[8:8 + sample_ct]
        phased, bits = block[8 + sample_ct], block[9 + sample_ct]
        if phased != 1:
            fail(varid + ": multiallelic block not in phased format")
        max_val = (1 << bits) - 1
        chrom = chroms[varid]
        is_x = (chrom in ("X", "23", "chrX")) and any_male
        is_y = chrom in ("Y", "24", "chrY")
        is_mt = chrom in ("MT", "M", "26", "chrM")
        tokens = vals[varid]
        expected_numers = []
        expected_ploidy = []
        for sidx in range(sample_ct):
            tok = tokens[sidx]
            missing = (tok == "M")
            u = d = None
            if not missing:
                parts = tok.split(";")
                u = [int(x) for x in parts[0].split(",")]
                if len(parts) > 1:
                    d = [int(x) for x in parts[1].split(",")]
            split = d is not None
            if is_x:
                ploidy = 1 if ((sexes[sidx] == "1") and not split) else 2
            elif is_y and any_female:
                ploidy = 0 if (missing and sexes[sidx] == "2") else (2 if split else 1)
            elif is_y or is_mt:
                ploidy = 2 if split else 1
            else:
                ploidy = 2
            expected_ploidy.append(ploidy + (128 if missing else 0))
            if missing:
                expected_numers.extend([0] * (ploidy * (allele_ct - 1)))
                continue
            if split:
                haps = [[u[i] + d[i] for i in range(allele_ct)],
                        [u[i] - d[i] for i in range(allele_ct)]]
                split_samples += 1
            else:
                haps = [u, u]
            if ploidy == 1:
                haploid_samples += 1
            for hap in haps[:ploidy]:
                if sum(hap) != 32768 or min(hap) < 0:
                    fail(varid + ": bad expected haplotype values")
                out_hap = [hap[out_to_in[i]] for i in range(allele_ct)]
                expected_numers.extend(round_hap(out_hap, max_val)[:-1])
            checked_samples += 1
        if list(ploidy_bytes) != expected_ploidy:
            for sidx in range(sample_ct):
                if ploidy_bytes[sidx] != expected_ploidy[sidx]:
                    fail("%s sample %d: ploidy/missingness byte %d, expected %d" % (varid, sidx, ploidy_bytes[sidx], expected_ploidy[sidx]))
        ploidies = [p & 63 for p in expected_ploidy]
        if (pmin, pmax) != (min(ploidies), max(ploidies)):
            fail("%s: ploidy range (%d, %d), expected (%d, %d)" % (varid, pmin, pmax, min(ploidies), max(ploidies)))
        prob_bytes = block[10 + sample_ct:]
        got = unpack_bits(prob_bytes, bits, len(expected_numers))
        if (len(prob_bytes) * 8 < len(expected_numers) * bits) or (len(prob_bytes) != (len(expected_numers) * bits + 7) // 8):
            fail(varid + ": probability block length mismatch")
        if got != expected_numers:
            for i in range(len(got)):
                if got[i] != expected_numers[i]:
                    fail("%s: probability #%d is %d, expected %d" % (varid, i, got[i], expected_numers[i]))
        checked_variants += 1
    if pos != len(raw):
        fail("trailing bytes")
    if ref_last_seen != set([1 if ref_mode == "ref-last" else 0]):
        fail("REF allele position doesn't match " + ref_mode)
    sample_fname = bgen_fname[:-len(".bgen")] + ".sample"
    sample_lines = [line.split() for line in open(sample_fname)][2:]
    if len(sample_lines) != sample_ct:
        fail(".sample row count mismatch")
    for sidx, f in enumerate(sample_lines):
        expected = missing_cts[sidx] / variant_ct
        if abs(float(f[2]) - expected) > 1e-5 * max(expected, 1e-3):
            fail("sample %d: .sample missing value %s, expected %g" % (sidx, f[2], expected))
    if not checked_variants:
        fail("no multiallelic variants checked")
    print("bgen_multiallelic_check: %s ok (%d multiallelic variants, %d sample calls, %d split, %d haploid)" % (bgen_fname, checked_variants, checked_samples, split_samples, haploid_samples))


main()

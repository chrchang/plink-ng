#!/usr/bin/env python3
"""Randomized oracle for multipass --pmerge-list.

The only hard requirement for a multipass merge is that its output matches the
corresponding single-pass merge byte for byte.  This builds random collections
of overlapping filesets (split by sample and by variant, with conflicting and
missing genotypes, phase, dosages, multiallelic variants, provisional REF
alleles, and conflicting QUAL/FILTER/INFO/CM values), merges each collection
once in a single pass and again with small forced pass sizes, and compares the
.pgen/.pvar/.psam outputs.

Usage: multipass_oracle.py <plink2 binary> <config count> <seed> [extra plink2
       args...]

Standard library only (the CI runners have no numpy).
"""

import os
import random
import shutil
import subprocess
import sys

CHROMS = ["1", "2"]
BASES = ["A", "C", "G", "T"]
QUAL_VALS = ["30", "30.0", "29.9999999", "30.00001", "7.25", "45"]
FILTER_VALS = ["PASS", "q10", "s50", "q10;s50", "s50;q10", "q10;q10", "PASS;q10"]
CM_VALS = ["0.5", "0.50", "0.5000000001", "1.25", "0"]
# key: (Number, Type)
INFO_KEYS = {
    "DP": ("1", "Integer"),
    "AF": ("A", "Float"),
    "RC": ("R", "Integer"),
    "FL": ("0", "Flag"),
    "TW": ("2", "Integer"),
    "XU": (".", "String"),
}


class Variant(object):
    def __init__(self, chrom, pos, vid, alleles, cm):
        self.chrom = chrom
        self.pos = pos
        self.vid = vid
        self.alleles = alleles
        self.cm = cm
        self.qual = None
        self.filt = None
        self.info = {}
        self.geno = {}
        self.dosage = {}


def make_universe(rng, sample_ids):
    variants = []
    for chrom in CHROMS:
        pos_ct = rng.randint(4, 9)
        positions = sorted(rng.sample(range(100, 400), pos_ct))
        for pos in positions:
            id_ct = rng.choice([1, 1, 1, 2])
            for id_idx in range(id_ct):
                # "rs9" vs. "rs10" orders differently under natural and ASCII
                # sorting.
                vid = "rs%d" % (pos * 10 + id_idx * 9 + (1 if id_idx else 0))
                ref = rng.choice(BASES)
                if rng.random() < 0.15:
                    ref += rng.choice(BASES)
                alt_pool = [b for b in BASES if b != ref] + [ref + "T", ref + "GA"]
                alt_ct = rng.choice([1, 1, 1, 2, 3])
                alts = rng.sample(alt_pool, alt_ct)
                variants.append(Variant(chrom, pos, vid, [ref] + alts, rng.choice(CM_VALS)))
    for v in variants:
        allele_ct = len(v.alleles)
        v.qual = rng.choice(QUAL_VALS)
        v.filt = rng.choice(FILTER_VALS)
        v.info = {
            "DP": rng.choice(["5", "6"]),
            "AF": [rng.choice(["0.1", "0.2", "0.10"]) for _ in range(allele_ct - 1)],
            "RC": [rng.choice(["1", "2"]) for _ in range(allele_ct)],
            "TW": rng.choice(["1,2", "2,1"]),
            "XU": rng.choice(["a", "b"]),
        }
        for sid in sample_ids:
            if allele_ct == 2:
                a1 = rng.randint(0, 1)
                a2 = rng.randint(0, 1)
            else:
                a1 = rng.randint(0, allele_ct - 1)
                a2 = rng.randint(0, allele_ct - 1)
            v.geno[sid] = (a1, a2, rng.random() < 0.6)
            if allele_ct == 2:
                # haplotype dosages
                v.dosage[sid] = (rng.choice(["0", "0.25", "1", "0.9"]), rng.choice(["0", "0.5", "1", "0.95"]))
    return variants


def pick(rng, truth, alternatives, p_truth, p_missing, missing="."):
    x = rng.random()
    if x < p_truth:
        return truth
    if x < p_truth + p_missing:
        return missing
    return rng.choice(alternatives)


class Fileset(object):
    pass


def make_fileset(rng, idx, variants, sample_ids, pr_mode, missing_ids, workdir, plink2, extra):
    fs = Fileset()
    fs.prefix = os.path.join(workdir, "fs%d" % idx)
    samples = rng.sample(sample_ids, rng.randint(1, len(sample_ids)))
    if rng.random() < 0.5:
        samples.sort(key=sample_ids.index)
    fs.samples = samples
    keep_frac = rng.choice([0.3, 0.6, 0.9, 1.0])
    recs = [v for v in variants if rng.random() < keep_frac]
    if not recs:
        recs = [rng.choice(variants)]
    dosage_mode = rng.choice(["none", "DS", "HDS"])
    if dosage_mode != "none":
        # VCF import doesn't support multiallelic dosages yet.
        recs = [v for v in recs if len(v.alleles) == 2]
        if not recs:
            dosage_mode = "none"
            recs = [rng.choice(variants)]
    has_qual = rng.random() < 0.7
    has_filter = rng.random() < 0.7
    has_info = rng.random() < 0.8
    has_cm = rng.random() < 0.5
    phased = rng.random() < 0.6
    if pr_mode == "info_pr":
        has_info = True
    p_missing = rng.choice([0.05, 0.2, 0.4])
    p_noise = rng.choice([0.0, 0.1, 0.3])
    pvar_lines = []
    vcf_lines = []
    info_pr_header = (pr_mode == "info_pr") or ((pr_mode == "all_pr") and has_info and (rng.random() < 0.5))
    conflicting_dp = rng.random() < 0.08
    for v in recs:
        allele_ct = len(v.alleles)
        is_pr = (pr_mode == "all_pr") or ((pr_mode == "info_pr") and (rng.random() < 0.35))
        # canonical allele indices, in this record's order
        if allele_ct == 2:
            order = [0, 1]
            if is_pr and rng.random() < 0.5:
                order = [1, 0]
        else:
            alt_idxs = list(range(1, allele_ct))
            rng.shuffle(alt_idxs)
            alt_idxs = alt_idxs[:rng.randint(1, len(alt_idxs))]
            order = [0] + alt_idxs
            if is_pr and rng.random() < 0.5:
                rng.shuffle(order)
        # A known-REF biallelic record may have a missing ALT allele, as when
        # it's derived from a .ped file.  Only REF calls are possible then.
        missing_alt = (allele_ct == 2) and (not is_pr) and (rng.random() < 0.08)
        if missing_alt:
            order = [0]
        pos_of = {canon: i for i, canon in enumerate(order)}
        rec_alleles = [v.alleles[c] for c in order]
        if missing_alt:
            rec_alleles.append(".")
        fmt_dosage = (dosage_mode != "none") and (allele_ct == 2) and (not missing_alt)
        fmt = "GT"
        if fmt_dosage:
            fmt += ":" + dosage_mode
        calls = []
        for sid in samples:
            a1, a2, ph = v.geno[sid]
            x = rng.random()
            if x < p_missing:
                gt = "./."
                call_alleles = None
            else:
                if (x < p_missing + p_noise) or (a1 not in pos_of) or (a2 not in pos_of):
                    a1 = rng.choice(order)
                    a2 = rng.choice(order)
                if missing_alt and (a1 or a2):
                    a1 = 0
                    a2 = 0
                sep = "|" if (phased and ph) else "/"
                gt = "%d%s%d" % (pos_of[a1], sep, pos_of[a2])
                call_alleles = (a1, a2)
            call = gt
            if fmt_dosage:
                if rng.random() < 0.5:
                    h1, h2 = v.dosage[sid]
                    if order == [1, 0]:
                        h1 = "%g" % (1 - float(h1))
                        h2 = "%g" % (1 - float(h2))
                    if dosage_mode == "HDS":
                        call += ":%s,%s" % (h1, h2)
                    else:
                        call += ":%g" % (float(h1) + float(h2))
                else:
                    call += ":."
            calls.append(call)
        vcf_lines.append("\t".join([v.chrom, str(v.pos), v.vid, rec_alleles[0], ",".join(rec_alleles[1:]), ".", ".", ".", fmt] + calls))
        # .pvar fields
        vid = v.vid
        if missing_ids and (rng.random() < 0.3):
            vid = "."
        fields = [v.chrom, str(v.pos), vid, rec_alleles[0], ",".join(rec_alleles[1:])]
        if has_qual:
            fields.append(pick(rng, v.qual, QUAL_VALS, 0.6, 0.25))
        if has_filter:
            fields.append(pick(rng, v.filt, FILTER_VALS, 0.6, 0.25))
        if has_info:
            kvs = []
            keys = list(INFO_KEYS.keys())
            rng.shuffle(keys)
            for key in keys:
                if rng.random() < 0.35:
                    continue
                if missing_alt and (key in ("AF", "RC")):
                    continue
                num = INFO_KEYS[key][0]
                if num == "0":
                    kvs.append(key)
                elif num == "A":
                    vals = [pick(rng, v.info["AF"][c - 1], ["0.3", "0.1"], 0.7, 0.2) for c in order[1:]] if 0 not in order[1:] else None
                    if vals is None:
                        # swapped PR record: the canonical REF is an ALT here
                        vals = [rng.choice(["0.3", "0.1", "."]) for _ in order[1:]]
                    kvs.append("AF=" + ",".join(vals))
                elif num == "R":
                    vals = [pick(rng, v.info["RC"][c], ["3", "1"], 0.7, 0.2) for c in order]
                    kvs.append("RC=" + ",".join(vals))
                elif num == "2":
                    kvs.append("TW=" + pick(rng, v.info["TW"], ["1,1", ".,."], 0.6, 0.0))
                elif key == "DP":
                    kvs.append("DP=" + pick(rng, v.info["DP"], ["7", "5"], 0.6, 0.2))
                else:
                    kvs.append("XU=" + pick(rng, v.info["XU"], ["c", "a"], 0.6, 0.2))
            if is_pr and info_pr_header and ((pr_mode == "info_pr") or (rng.random() < 0.5)):
                kvs.insert(rng.randint(0, len(kvs)), "PR")
            fields.append(";".join(kvs) if kvs else ".")
        if has_cm:
            fields.append(pick(rng, v.cm, CM_VALS, 0.6, 0.0))
        pvar_lines.append("\t".join(fields))
    vcf_path = fs.prefix + ".vcf"
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for chrom in CHROMS:
            f.write("##contig=<ID=%s,length=1000>\n" % chrom)
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">\n')
        f.write('##FORMAT=<ID=DS,Number=A,Type=Float,Description="ds">\n')
        f.write('##FORMAT=<ID=HDS,Number=.,Type=Float,Description="hds">\n')
        f.write("#" + "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + ["s%d" % s for s in samples]) + "\n")
        for line in vcf_lines:
            f.write(line + "\n")
    cmd = [plink2] + extra + ["--vcf", vcf_path]
    if dosage_mode != "none":
        cmd.append("dosage=" + dosage_mode)
    cmd += ["--double-id", "--make-pgen", "--out", fs.prefix]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # Replace the imported .pvar with our own.  Genotype indexes are unchanged
    # since every record keeps its allele count.
    with open(fs.prefix + ".pvar", "w") as f:
        for key, (num, typ) in INFO_KEYS.items():
            if conflicting_dp and (key == "DP"):
                typ = "Float"
            f.write('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % (key, num, typ, key))
        if info_pr_header:
            f.write('##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">\n')
        f.write('##FILTER=<ID=q10,Description="q10">\n##FILTER=<ID=s50,Description="s50">\n')
        cols = ["#CHROM", "POS", "ID", "REF", "ALT"]
        if has_qual:
            cols.append("QUAL")
        if has_filter:
            cols.append("FILTER")
        if has_info:
            cols.append("INFO")
        if has_cm:
            cols.append("CM")
        f.write("\t".join(cols) + "\n")
        for line in pvar_lines:
            f.write(line + "\n")
    if pr_mode == "all_pr":
        # Flip the .pgen's provisional-REF storage mode from "always trusted"
        # to "always untrusted".  (Both modes store no per-variant flags, so
        # the header layout is unchanged.)
        with open(fs.prefix + ".pgen", "r+b") as f:
            f.seek(11)
            ctrl = f.read(1)[0]
            assert (ctrl >> 6) == 1, ctrl
            f.seek(11)
            f.write(bytes([(ctrl & 0x3f) | 0x80]))
    return fs


def read_bytes(path):
    with open(path, "rb") as f:
        return f.read()


def run_merge(plink2, extra, list_path, out, flags, pass_size):
    cmd = [plink2] + extra + ["--pmerge-list", list_path] + flags + ["--out", out]
    if pass_size:
        cmd += ["--pmerge-pass-size", str(pass_size)]
    return subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode


def main():
    plink2 = sys.argv[1]
    config_ct = int(sys.argv[2])
    seed = int(sys.argv[3])
    extra = sys.argv[4:]
    rng = random.Random(seed)
    compared_ct = 0
    error_ct = 0
    pass_run_ct = 0
    multi_run_ct = 0
    error_msgs = {}
    concat_ct = 0
    for config_idx in range(config_ct):
        workdir = "tmp_cfg%d" % config_idx
        if os.path.exists(workdir):
            shutil.rmtree(workdir)
        os.mkdir(workdir)
        sample_ids = list(range(1, rng.randint(2, 14)))
        variants = make_universe(rng, sample_ids)
        # Mostly a few filesets (forced small passes), sometimes more than 20
        # so that the default pass size kicks in by itself.
        fileset_ct = rng.choice([2, 3, 4, 5, 6, 7, 9, 12]) if (config_idx % 5) else rng.randint(21, 45)
        pr_modes = ["trusted", "trusted", "info_pr", "all_pr"]
        if rng.random() < 0.3:
            pr_modes = ["trusted"]
        missing_ids = rng.random() < 0.15
        filesets = [make_fileset(rng, i, variants, sample_ids, rng.choice(pr_modes), missing_ids, workdir, plink2, extra) for i in range(fileset_ct)]
        list_path = os.path.join(workdir, "list.txt")
        with open(list_path, "w") as f:
            for fs in filesets:
                f.write(os.path.basename(fs.prefix) + "\n")
        flags = ["--pmerge-list-dir", workdir]
        flags += ["--merge-mode", rng.choice(["nm-match", "nm-first", "first"])]
        flags += ["--merge-qual-mode", rng.choice(["nm-match", "nm-first", "first", "min"])]
        flags += ["--merge-filter-mode", rng.choice(["nm-match", "nm-first", "first", "np-union"])]
        flags += ["--merge-info-mode", rng.choice(["nm-match", "nm-first", "first"])]
        flags += ["--merge-cm-mode", rng.choice(["nm-match", "nm-first", "first"])]
        if rng.random() < 0.3:
            flags += ["--merge-info-sort", rng.choice(["ascii", "natural"])]
        if rng.random() < 0.2:
            flags += ["--sort-vars", "ascii"]
        if missing_ids:
            flags += ["--set-missing-var-ids", "@:#:$r:$a"]
        if rng.random() < 0.15:
            flags.append("--merge-ignore-phase")
        elif rng.random() < 0.15:
            flags.append("--merge-ignore-dosage")
        if (rng.random() < 0.15) and set.intersection(*[set(fs.samples) for fs in filesets]):
            flags.append("--sample-inner-join")
        if rng.random() < 0.15:
            flags.append("--pmerge-output-vzs")
        single_out = os.path.join(workdir, "single")
        # The single-pass reference: every fileset merged at once.
        single_rc = run_merge(plink2, extra, list_path, single_out, flags, 82)
        pvar_ext = ".pvar.zst" if "--pmerge-output-vzs" in flags else ".pvar"
        if single_rc == 0:
            ref = [read_bytes(single_out + ext) for ext in (".pgen", pvar_ext, ".psam")]
        else:
            with open(single_out + ".log") as f:
                errors = [line[line.index("Error"):].strip() for line in f if "Error" in line]
            error_msgs[errors[0] if errors else "(no message)"] = error_msgs.get(errors[0] if errors else "(no message)", 0) + 1
        pass_sizes = sorted(set([2, rng.randint(2, 4), rng.randint(2, max(2, fileset_ct - 1))]))
        if fileset_ct > 20:
            pass_sizes.append(0)
        with open(single_out + ".log") as f:
            if "Concatenation job detected" in f.read():
                # Rare: the filesets happen not to overlap, so there's only
                # ever one pass.
                concat_ct += 1
                shutil.rmtree(workdir)
                continue
        for pass_size in pass_sizes:
            if pass_size >= fileset_ct:
                continue
            multi_out = os.path.join(workdir, "multi%d" % pass_size)
            rc = run_merge(plink2, extra, list_path, multi_out, flags, pass_size)
            multi_run_ct += 1
            if (rc == 0) != (single_rc == 0):
                sys.stderr.write("config %d (%s): single-pass exit code %d, pass size %d exit code %d\n" % (config_idx, " ".join(flags), single_rc, pass_size, rc))
                sys.exit(1)
            if rc:
                error_ct += 1
                continue
            for ext, ref_bytes in zip((".pgen", pvar_ext, ".psam"), ref):
                if read_bytes(multi_out + ext) != ref_bytes:
                    sys.stderr.write("config %d (%s): pass size %d %s differs from single-pass output\n" % (config_idx, " ".join(flags), pass_size, ext))
                    sys.exit(1)
            with open(multi_out + ".log") as f:
                if " pass 1: " not in f.read():
                    sys.stderr.write("config %d: pass size %d didn't produce a multipass merge\n" % (config_idx, pass_size))
                    sys.exit(1)
            leftovers = [fn for fn in os.listdir(workdir) if "-merge-tmp" in fn]
            if leftovers:
                sys.stderr.write("config %d: temporary files left behind: %s\n" % (config_idx, " ".join(leftovers)))
                sys.exit(1)
            compared_ct += 1
        pass_run_ct += 1
        shutil.rmtree(workdir)
    print("%d configurations (plus %d concatenation jobs), %d multipass merges compared byte-for-byte, %d merges failed in both modes." % (pass_run_ct, concat_ct, compared_ct, error_ct))
    for msg, ct in sorted(error_msgs.items()):
        print("  single-pass error (%d configurations): %s" % (ct, msg))


if __name__ == "__main__":
    main()

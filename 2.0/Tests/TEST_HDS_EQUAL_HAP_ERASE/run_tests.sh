#!/bin/bash

# --dosage-erase-threshold on phased dosages whose haplotype dosages are equal
# (h|h, or a phased haploid h) must use the per-haplotype rule (erase only if
# every haplotype is within threshold/2 of an integer), in both the scan and
# the conversion pass, for VCF HDS and .bgen import alike.
#
# .bgen import used to apply the unphased rule here (with
# --dosage-erase-threshold 0.05, 0.498|0.498 = dosage 0.996 is within 0.05 of
# 1, so it was erased to a het hardcall), while VCF/BCF HDS conversion applied
# the per-haplotype rule and kept it.  The .bgen scan pass also has to agree
# with conversion: if it treats a record as dosage-free and nothing else in the
# file has dosages, the dosage track is never written.  So each kind of record
# is imported both alone and next to a record that forces the dosage path, and
# the results must match.
#
# Records (samples m = male, f = female):
#   v1 (chr1):  0.498|0.498, 0|0              fixed ploidy 2
#   v2 (chrX):  0.498 (haploid), 0|0          variable ploidy
#   v4 (chrMT): 0.498, 0 (both haploid)       fixed ploidy 1
#   v3 (chrX):  forces the dosage path (0.3|0.7)

set -exo pipefail

printf '#IID\tSEX\nm\t1\nf\t2\n' > tmp.psam

# $1: output name, remaining args: record IDs
write_vcf() {
    local out=$1 rec
    shift
    {
        printf '##fileformat=VCFv4.3\n##contig=<ID=1,length=1000000>\n##contig=<ID=X,length=156040895>\n##contig=<ID=MT,length=16569>\n'
        printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        printf '##FORMAT=<ID=HDS,Number=.,Type=Float,Description="Haplotype dosages">\n'
        printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tm\tf\n'
        for rec in "$@"; do
            case $rec in
                v1) printf '1\t100\tv1\tA\tC\t.\t.\t.\tGT:HDS\t0|0:0.498,0.498\t0|0:0,0\n' ;;
                v2) printf 'X\t5000000\tv2\tA\tC\t.\t.\t.\tGT:HDS\t0:0.498\t0|0:0,0\n' ;;
                v3) printf 'X\t6000000\tv3\tA\tC\t.\t.\t.\tGT:HDS\t0:0\t0|1:0.3,0.7\n' ;;
                v4) printf 'MT\t300\tv4\tA\tC\t.\t.\t.\tGT:HDS\t0:0.498\t0:0\n' ;;
            esac
        done
    } > $out
}

# BGEN versions at 16 bits (stored value = P(haplotype carries REF);
# 0.502 * 65535 rounds to 32899, i.e. ALT dosage 0.49799).
write_bgen() {
    python3 - "$@" << 'PY'
import os, sys
sys.path.insert(0, os.path.join("..", "utils"))
import bgen_writer
d = 65535
h = round(0.502 * d)
def var(chrom, pos, vid, ploidy, numers):
    return {"chrom": chrom, "pos": pos, "varid": vid, "rsid": vid,
            "alleles": ["A", "C"], "phased": True, "ploidy": ploidy,
            "missing": [False] * len(ploidy), "numers": numers}
recs = {
    "v1": var("1", 100, "v1", [2, 2], [[h, h], [d, d]]),
    "v2": var("X", 5000000, "v2", [1, 2], [[h], [d, d]]),
    "v3": var("X", 6000000, "v3", [1, 2], [[d], [round(0.7 * d), round(0.3 * d)]]),
    "v4": var("MT", 300, "v4", [1, 1], [[h], [d]]),
}
bgen_writer.write_bgen(sys.argv[1], ["m", "f"], [recs[v] for v in sys.argv[2:]],
                       2, 16, "zlib", embed_ids=False)
PY
}
printf 'ID_1 ID_2 missing sex\n0 0 0 D\nm m 0 1\nf f 0 2\n' > tmp.sample

# Prints the given record's calls (GT, DS, HDS per sample).  Equal haplotypes
# carry no phase, and '|' vs '/' on the other calls depends on whether
# anything else in the file is phased, so separators are normalized.
summarize() {
    $1/plink2 $2 $3 --pfile $4 --export vcf vcf-dosage=HDS-force --out $4_e > /dev/null
    grep -v '^#' $4_e.vcf | awk -v v=$5 '$3 == v' | cut -f 3,10- | tr '|' '/'
}

for v in v1 v2 v4; do
    write_vcf tmp_${v}_alone.vcf $v
    write_vcf tmp_${v}_forced.vcf $v v3
    write_bgen tmp_${v}_alone.bgen $v
    write_bgen tmp_${v}_forced.bgen $v v3
    for fmt in vcf bgen; do
        for f in alone forced; do
            if [ $fmt = vcf ]; then
                src="--vcf tmp_${v}_$f.vcf dosage=HDS --psam tmp.psam"
            else
                src="--bgen tmp_${v}_$f.bgen ref-first --sample tmp.sample"
            fi
            $1/plink2 $2 $3 $src --dosage-erase-threshold 0.05 --make-pgen --out tmp_${v}_${fmt}_$f
        done
        summarize $1 "$2" "$3" tmp_${v}_${fmt}_alone $v > tmp_${v}_${fmt}_alone.txt
        summarize $1 "$2" "$3" tmp_${v}_${fmt}_forced $v > tmp_${v}_${fmt}_forced.txt
        diff tmp_${v}_${fmt}_alone.txt tmp_${v}_${fmt}_forced.txt
        # The dosage must survive: each haplotype at ~0.498 is far from an
        # integer.
        grep -q $'^'$v$'\t[^\t]*:0.996:0.498' tmp_${v}_${fmt}_alone.txt
    done
    diff tmp_${v}_vcf_alone.txt tmp_${v}_bgen_alone.txt
done

#!/bin/bash

# --set-invalid-haploid-missing / --set-mixed-mt-missing with dosages.
#
# After clearing the first affected call's dosage, the remaining dosages were
# not shifted down, so later samples got the wrong dosage values.  The result
# must equal importing the same data with those calls already missing.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
import random
random.seed(41)
n = 40
ids = ["s%d" % i for i in range(n)]
sexes = [1 if i % 4 else 2 for i in range(n)]
with open("tmp_data.psam", "w") as f:
    f.write("#IID\tSEX\n")
    for i in range(n):
        f.write("%s\t%d\n" % (ids[i], sexes[i]))
def write(fn, clean):
    with open(fn, "w") as f:
        f.write("##fileformat=VCFv4.2\n##contig=<ID=Y,length=50000000>\n##contig=<ID=MT,length=16569>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=DS,Number=A,Type=Float,Description="Dosage">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
        random.seed(42)
        for chrom, pos0 in (("Y", 5000000), ("MT", 1000)):
            for v in range(4):
                cells = []
                for i in range(n):
                    r = random.random()
                    female_y = (chrom == "Y" and sexes[i] == 2)
                    if r < 0.3:
                        cell = "0/1:%.3f" % (0.9 + 0.2 * random.random())
                        bad = True
                    else:
                        a = int(random.random() < 0.5)
                        cell = "%d:%.3f" % (a, a * 0.9 + 0.05 * random.random())
                        bad = False
                    if clean and (bad or female_y):
                        cell = "."
                    cells.append(cell)
                f.write("%s\t%d\t%s%d\tA\tG\t.\t.\t.\tGT:DS\t%s\n" % (chrom, pos0 + 100 * v, chrom, v, "\t".join(cells)))
write("tmp_data.vcf", False)
write("tmp_clean.vcf", True)
PY
$plink2 --vcf tmp_data.vcf dosage=DS --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null
$plink2 --pfile tmp_data --set-invalid-haploid-missing --set-mixed-mt-missing --make-pgen --out tmp_set > /dev/null
$plink2 --vcf tmp_clean.vcf dosage=DS --psam tmp_data.psam --make-pgen --out tmp_clean > /dev/null
$plink2 --pfile tmp_set --export vcf vcf-dosage=DS --out tmp_set > /dev/null
$plink2 --pfile tmp_clean --export vcf vcf-dosage=DS --out tmp_clean > /dev/null
diff <(grep -v '^##' tmp_set.vcf) <(grep -v '^##' tmp_clean.vcf)

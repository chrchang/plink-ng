#!/bin/bash

# --score 'center' and 'variance-standardize' with missing genotypes.
#
# Both modifiers shift every genotype so that its expected value is zero.  A
# missing genotype is mean-imputed to the expected value, so after the shift it
# must contribute nothing at all: the score sums with and without
# 'no-mean-imputation' have to agree.  The imputed value used to be left
# unshifted, adding 2 * freq * weight (divided by the standard deviation under
# 'variance-standardize') for every missing call.
#
# Two datasets exercise the three code paths that fill in missing genotypes:
# rare variants stored as sparse difflists and common variants read as plain
# hardcalls in the first, variants with dosages in the second.  The second
# also has chrY and chrM variants, which are haploid.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# $1: output prefix.  $2: 1 to add DS dosages.
make_vcf() {
    awk -v with_ds=$2 'BEGIN {
      srand(17 + with_ds); OFS = "\t"; n = 203
      print "##fileformat=VCFv4.2"
      print "##contig=<ID=1,length=100000000>"
      print "##contig=<ID=Y,length=60000000>"
      print "##contig=<ID=MT,length=16569>"
      print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
      if (with_ds) { print "##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"Alternate allele dosage\">" }
      hdr = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
      for (s = 1; s <= n; ++s) { hdr = hdr "\ts" s }
      print hdr
      nvar = 40
      for (v = 1; v <= nvar; ++v) {
        chrom = "1"; pos = v * 1000
        if (with_ds && (v > 30)) { chrom = (v > 35)? "MT" : "Y"; pos = 3000000 + v * 10 }
        haploid = (chrom != "1")
        # rare alleles for the first dozen variants, common ones afterwards
        f = (v <= 12)? 0.004 + 0.002 * (v % 3) : 0.1 + 0.8 * rand()
        line = chrom OFS pos OFS "v" v OFS "A" OFS "G" OFS "." OFS "." OFS "." OFS (with_ds? "GT:DS" : "GT")
        for (s = 1; s <= n; ++s) {
          # a missing call for every variant in about 1 sample in 12
          if (rand() < 0.08) {
            g = haploid? "." : "./."
            if (with_ds) { g = g ":." }
          } else {
            a = (rand() < f); b = (rand() < f)
            if (haploid) { g = a; d = a } else { g = a "/" b; d = a + b }
            if (with_ds) {
              if (rand() < 0.3) { d = d + (rand() - 0.5) * 0.6; if (d < 0) { d = 0 }; if (d > 2 - haploid) { d = 2 - haploid } }
              g = g ":" sprintf("%.3f", d)
            }
          }
          line = line OFS g
        }
        print line
      }
    }' > $1.vcf
}

# The score file: every variant, with two weight columns.  The rare variants
# are scored on ALT; the rest alternate between ALT and REF.  (A rare variant
# can come out monomorphic, and 'variance-standardize' rejects scoring the
# REF allele of a variant with no ALT call.)
make_scores() {
    grep -v '^#' $1.pvar | awk 'BEGIN { srand(5); OFS = "\t" } { print $3, ((NR <= 12) || (NR % 2))? $5 : $4, sprintf("%.4f", rand() * 2 - 1), sprintf("%.4f", rand() * 4 - 2) }' > $1.score
}

# $1: fileset prefix.  $2: the modifier.  $3: any extra flags.
check_mode() {
    $plink2 --pfile $1 --score $1.score 1 2 $2 cols=+scoresums --score-col-nums 3,4 $3 --out tmp_imp > /dev/null
    $plink2 --pfile $1 --score $1.score 1 2 $2 no-mean-imputation cols=+scoresums --score-col-nums 3,4 $3 --out tmp_noimp > /dev/null
    # SCORE1_SUM and SCORE2_SUM are the last two columns.  Values are printed
    # to six significant digits, and the two runs add their terms in different
    # orders, so allow a small absolute slack.
    awk 'function abs(x) { return (x < 0)? -x : x }
         FNR == 1 { next }
         NR == FNR { s1[$1] = $(NF - 1); s2[$1] = $NF; next }
         {
           if (!($1 in s1)) { print "sample " $1 " missing"; exit 1 }
           if ((abs(s1[$1] - $(NF - 1)) > 1e-4) || (abs(s2[$1] - $NF) > 1e-4)) { print "sample " $1 ": " s1[$1] " " s2[$1] " vs " $(NF - 1) " " $NF; exit 1 }
           ++n
         }
         END { if (n != 203) { print "expected 203 samples, got " n + 0; exit 1 } }' tmp_imp.sscore tmp_noimp.sscore
}

# 1. Hardcalls only; the rare variants are stored as difflists.
make_vcf tmp_hc 0
$plink2 --vcf tmp_hc.vcf --make-pgen --out tmp_hc > /dev/null
make_scores tmp_hc
check_mode tmp_hc center
check_mode tmp_hc variance-standardize

# 2. The same with some variants given dosages, plus chrY and chrM.  Every
#    sample is male, so chrY is scored for all of them.
make_vcf tmp_ds 1
awk 'BEGIN { print "#IID\tSEX"; for (s = 1; s <= 203; ++s) { print "s" s "\t1" } }' > tmp_ds.sex
$plink2 --vcf tmp_ds.vcf dosage=DS --update-sex tmp_ds.sex --make-pgen --out tmp_ds > /dev/null
make_scores tmp_ds
check_mode tmp_ds center
# ('variance-standardize' is not allowed on chrM.)
check_mode tmp_ds variance-standardize "--not-chr MT"

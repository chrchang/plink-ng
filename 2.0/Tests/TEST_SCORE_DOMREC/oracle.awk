# Expected --score SCORE1_SUM for each sample, computed from a VCF.
#
# Usage: awk -v model={additive|dominant|recessive} -v meanimpute={0|1}
#            -f oracle.awk <.afreq> <score file> <VCF>
#
# The dosage of the scored allele is taken from the DS field when the VCF has
# one, and from GT otherwise.  'dominant' codes it as min(dosage, 1) and
# 'recessive' as max(dosage - 1, 0).  A missing genotype contributes the
# expected coded value under Hardy-Weinberg equilibrium: 2f, 1 - (1-f)^2 or
# f^2, where f is the frequency of the scored allele in the .afreq file.

function coded(dosage) {
  if (model == "dominant") {
    return (dosage > 1)? 1 : dosage
  }
  if (model == "recessive") {
    return (dosage > 1)? dosage - 1 : 0
  }
  return dosage
}

function expected(f) {
  if (model == "dominant") {
    return 1 - (1 - f) * (1 - f)
  }
  if (model == "recessive") {
    return f * f
  }
  return 2 * f
}

BEGIN { OFS = "\t" }

FILENAME == ARGV[1] {
  if ($1 !~ /^#/) {
    ref[$2] = $3
    alt_freq[$2] = $5
  }
  next
}

FILENAME == ARGV[2] {
  score_allele[$1] = $2
  weight[$1] = $3
  next
}

/^##/ { next }

/^#CHROM/ {
  for (s = 10; s <= NF; ++s) {
    iid[s] = $s
  }
  last_col = NF
  next
}

{
  id = $3
  if (!(id in weight)) {
    next
  }
  is_ref = (score_allele[id] == ref[id])
  f = is_ref? 1 - alt_freq[id] : alt_freq[id]
  ds_idx = 0
  n_fmt = split($9, fmt, ":")
  for (k = 1; k <= n_fmt; ++k) {
    if (fmt[k] == "DS") {
      ds_idx = k
    }
  }
  for (s = 10; s <= NF; ++s) {
    split($s, field, ":")
    missing = 0
    if (ds_idx) {
      if (field[ds_idx] == ".") {
        missing = 1
      } else {
        alt_dosage = field[ds_idx] + 0
      }
    } else {
      gt = field[1]
      if (gt ~ /\./) {
        missing = 1
      } else {
        alt_dosage = gsub(/1/, "1", gt)
      }
    }
    if (missing) {
      if (meanimpute) {
        sum[s] += weight[id] * expected(f)
      }
    } else {
      sum[s] += weight[id] * coded(is_ref? 2 - alt_dosage : alt_dosage)
    }
  }
}

END {
  for (s = 10; s <= last_col; ++s) {
    printf("%s\t%.10g\n", iid[s], sum[s])
  }
}

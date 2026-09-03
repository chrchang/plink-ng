#!/bin/bash

# PLINK 1.9 --bmerge keys the merge on variant ID, so several records sharing an
# ID describe a single merged variant.  The log must therefore describe a
# fileset the same way whether it is the base or the second fileset, and the
# duplicate-ID warning must name whichever fileset has the duplicates.
# (https://github.com/chrchang/plink-ng/issues/140)

set -exo pipefail

# 4 records, 4 distinct IDs.
cat > plus.map << 'MAPEOF'
1	v1	0	100
1	v2	0	200
1	v3	0	300
1	v4	0	400
MAPEOF
cat > plus.ped << 'PEDEOF'
F1 S1 0 0 1 1 A G A A C C G G
F1 S2 0 0 1 1 A A A G C T G T
F1 S3 0 0 1 1 G G G G T T T T
PEDEOF

# First 3 variants of plus, with v2 present twice: 4 records, 3 distinct IDs.
cat > normed.map << 'MAPEOF'
1	v1	0	100
1	v2	0	200
1	v2	0	200
1	v3	0	300
MAPEOF
cat > normed.ped << 'PEDEOF'
F1 S1 0 0 1 1 A G A A A A C T
F1 S2 0 0 1 1 A A A G A G C T
F1 S3 0 0 1 1 G G G G G G T T
PEDEOF

# Same three variants, without the duplicate record.
sed -n '1,2p;4p' normed.map > single.map
cut -d ' ' -f 1-10,13-14 normed.ped > single.ped

plink --file plus --make-bed --out plus
plink --file normed --make-bed --out normed
plink --file single --make-bed --out single

# 1. Duplicate IDs in the second fileset: the warning names that fileset, and
#    its variant count is the number of distinct IDs, not the record count.
plink --bfile plus --bmerge normed.bed normed.bim normed.fam --merge-mode 6 --out base_plus
grep -q "normed.bim contains duplicate variant ID(s)" base_plus.log
grep -q "^3 markers to be merged from normed.bim.$" base_plus.log
grep -q "^Of these, 0 are new, while 3 are present in the base dataset.$" base_plus.log

# 2. Same pair the other way around: the same fileset reports the same count.
plink --bfile normed --bmerge plus.bed plus.bim plus.fam --merge-mode 6 --out base_normed
grep -q "normed.bim contains duplicate variant ID(s)" base_normed.log
grep -q "^3 markers loaded from normed.bim.$" base_normed.log
grep -q "^4 markers to be merged from plus.bim.$" base_normed.log
grep -q "^Of these, 1 is new, while 3 are present in the base dataset.$" base_normed.log

# 3. Either way around, the merged fileset has the 4 distinct variants.
plink --bfile plus --bmerge normed.bed normed.bim normed.fam --make-bed --out merged_plus
plink --bfile normed --bmerge plus.bed plus.bim plus.fam --make-bed --out merged_normed
[[ "$(wc -l < merged_plus.bim)" -eq 4 ]]
[[ "$(wc -l < merged_normed.bim)" -eq 4 ]]

# 4. No duplicate IDs anywhere: no warning.
plink --bfile plus --bmerge single.bed single.bim single.fam --merge-mode 6 --out no_dups
! grep -q "contains duplicate variant ID" no_dups.log
grep -q "^3 markers to be merged from single.bim.$" no_dups.log

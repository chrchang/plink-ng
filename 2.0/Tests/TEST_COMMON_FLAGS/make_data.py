#!/usr/bin/env python3
"""Writes tmp_data.ped/.map (20 trios, 60 unrelated, chr1) and tmp_sex.ped/.map (chrX), plus what is known about them by construction."""
import random

rng = random.Random(20261103)
VARIANT_CT = 300
samples = []  # (fid, iid, father, mother, sex)
for i in range(20):
    samples += [(f"t{i}", "pa", "0", "0", 1), (f"t{i}", "ma", "0", "0", 2), (f"t{i}", "kid", "pa", "ma", 1 + i % 2)]
samples += [(f"u{i}", "s", "0", "0", 1 + i % 2) for i in range(60)]
index = {(s[0], s[1]): i for i, s in enumerate(samples)}

columns = []  # one genotype list per variant; a genotype is a (0/1, 0/1) ALT-allele pair
map_lines = []
mendel_truth = []
for v in range(VARIANT_CT):
    name = f"v{v}"
    map_lines.append(f"1 {name} 0 {(v + 1) * 1000}\n")
    if v % 40 == 7:
        # Everybody heterozygous: fails an HWE test, and every parent-kid pair is consistent.
        columns.append([(0, 1)] * len(samples))
        continue
    if v % 10 == 1:
        # Nearly a copy of the previous variant, so that the pair is in strong LD.
        col = list(columns[-1])
        for i in rng.sample(range(len(samples)), 2):
            col[i] = (0, 0) if col[i] == (0, 1) else (0, 1)
    else:
        af = rng.uniform(0.2, 0.5)
        col = [(int(rng.random() < af), int(rng.random() < af)) for _ in samples]
    for s, (fid, iid, father, mother, _) in enumerate(samples):
        if father != "0":
            pa = col[index[(fid, father)]]
            ma = col[index[(fid, mother)]]
            col[s] = (rng.choice(pa), rng.choice(ma))
    # A kid who is hom-ALT with a hom-REF father, or the reverse, is always an error.
    if v % 25 == 3:
        t = v % 20
        col[3 * t] = (0, 0)
        col[3 * t + 1] = (0, 0)
        col[3 * t + 2] = (1, 1)
        mendel_truth.append(f"t{t} v{v}\n")
    columns.append(col)

with open("tmp_data.map", "w") as f:
    f.writelines(map_lines)
with open("tmp_data.ped", "w") as f:
    for s, (fid, iid, father, mother, sex) in enumerate(samples):
        alleles = " ".join("AC"[g] for col in columns for g in col[s])
        f.write(f"{fid} {iid} {father} {mother} {sex} -9 {alleles}\n")
with open("tmp_mendel_truth.txt", "w") as f:
    f.writelines(sorted(mendel_truth))
with open("tmp_founders.txt", "w") as f:
    f.writelines(f"{s[0]} {s[1]}\n" for s in samples if s[2] == "0")
with open("tmp_ld_pairs.txt", "w") as f:
    f.writelines(f"v{v - 1} v{v}\n" for v in range(VARIANT_CT) if v % 10 == 1 and (v - 1) % 40 != 7)

# chrX: males hemizygous, females often heterozygous.  The first two samples
# carry the wrong sex in the input, and the third has none.
true_sex = [1 + i % 2 for i in range(60)]
given_sex = list(true_sex)
given_sex[0], given_sex[1], given_sex[2] = 2, 1, 0
xcols = []
for v in range(200):
    af = rng.uniform(0.2, 0.5)
    col = []
    for i in range(60):
        a = int(rng.random() < af)
        col.append((a, a) if true_sex[i] == 1 else (a, int(rng.random() < af)))
    xcols.append(col)
with open("tmp_sex.map", "w") as f:
    f.writelines(f"X x{v} 0 {(v + 1) * 5000000}\n" for v in range(200))
with open("tmp_sex.ped", "w") as f:
    for i in range(60):
        alleles = " ".join("AC"[g] for col in xcols for g in col[i])
        f.write(f"x{i} x{i} 0 0 {given_sex[i]} -9 {alleles}\n")
with open("tmp_sex_truth.txt", "w") as f:
    f.writelines(f"x{i} x{i} {true_sex[i]}\n" for i in range(60))

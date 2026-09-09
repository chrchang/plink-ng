#!/usr/bin/env python3
"""Independent reference implementation of --dfam, for testing."""
import sys
from collections import OrderedDict

MENDEL_TABLE = [
  0, 0, 0x6000101, 0,
  0, 0, 0x6000101, 0,
  0x7010001, 0x7010001, 0x8000001, 0x7010001,
  0, 0, 0x6000101, 0,
  0x2010101, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0x1010101, 0,
  0, 0, 0, 0,
  0x5000001, 0x4010001, 0x4010001, 0x4010001,
  0x3000101, 0, 0, 0,
  0x3000101, 0, 0, 0,
  0x3000101, 0, 0, 0]

# (dad_geno, mom_geno) -> total ALT count when informative, else 0.
PARENTAL = {(0, 1): 1, (1, 0): 1, (1, 1): 2, (1, 2): 3, (2, 1): 3}


def sibship_calc(case_ct, case_hom, case_het, ctrl_ct, ctrl_hom, ctrl_het,
                 plink19):
    """Returns (total_ct increment, numer, denom, expected)."""
    if not ctrl_ct:
        return (0, 0.0, 0.0, 0.0)
    hom = case_hom + ctrl_hom
    het = case_het + ctrl_het
    total = case_ct + ctrl_ct
    case_alt = 2 * case_hom + case_het
    if ((not hom) and (not het)) or het == total or hom == total:
        return (case_alt, 0.0, 0.0, float(case_alt))
    p = case_ct / total
    k = p * ctrl_ct / (total * (total - 1))
    var_hom = k * hom * (total - hom)
    var_het = k * het * (total - het)
    neg_covar = k * hom * het
    expected = 2 * p * hom + p * het
    # PLINK 1.x adds the covariance magnitude instead of subtracting it.
    var = (4 * (var_hom + neg_covar) + var_het if plink19 else
           4 * (var_hom - neg_covar) + var_het)
    return (case_alt, case_alt - expected, var, expected)


def main():
    psam_path, vcf_path = sys.argv[1], sys.argv[2]
    plink19 = len(sys.argv) > 3 and sys.argv[3] == 'plink19'
    # PLINK 1.x counts the minor allele rather than ALT, and its statistic is
    # not invariant under that choice, so its own reported A1 is used to pick
    # the orientation when reproducing it.
    counted_allele = {}
    if len(sys.argv) > 4:
        with open(sys.argv[4]) as f:
            for line in f:
                g = line.split()
                counted_allele[g[0]] = g[1]
    samples = OrderedDict()
    with open(psam_path) as f:
        header = f.readline().split()
        cols = {name.lstrip('#'): i for i, name in enumerate(header)}
        for line in f:
            g = line.split()
            iid = g[cols['IID']]
            samples[iid] = {'fid': g[cols['FID']], 'pat': g[cols['PAT']],
                            'mat': g[cols['MAT']], 'pheno': g[cols['PHENO1']]}
    # Only samples with a case/control phenotype take part.
    kept = [iid for iid in samples if samples[iid]['pheno'] in ('1', '2')]
    kept_set = set(kept)
    case = {iid: samples[iid]['pheno'] == '2' for iid in kept}
    founder = {iid: samples[iid]['pat'] == '0' and samples[iid]['mat'] == '0'
               for iid in kept}

    # Families: children whose parents are both present.  Keyed by parent pair,
    # in first-child order.
    families = OrderedDict()
    for iid in kept:
        pat, mat = samples[iid]['pat'], samples[iid]['mat']
        if pat in kept_set and mat in kept_set:
            families.setdefault((pat, mat), []).append(iid)
    in_family = set()
    for (pat, mat), children in families.items():
        in_family.add(pat)
        in_family.add(mat)
        in_family.update(children)

    all_case_groups = []
    mixed_groups = []
    for (pat, mat), children in families.items():
        case_ct = sum(case[c] for c in children)
        if case_ct == len(children):
            all_case_groups.append((pat, mat, children))
        elif case_ct:
            mixed_groups.append((pat, mat, children))

    # Sibships: non-founders who aren't a child in one of those families.
    trio_children = set()
    for children in families.values():
        trio_children.update(children)
    sib_keys = OrderedDict()
    for iid in kept:
        if founder[iid] or iid in trio_children:
            continue
        key = (samples[iid]['fid'], samples[iid]['pat'], samples[iid]['mat'])
        sib_keys.setdefault(key, []).append(iid)
    in_sibship = set()
    sibship_groups = []
    for key in sorted(sib_keys):
        members = sib_keys[key]
        if len(members) < 2:
            continue
        in_sibship.update(members)
        case_ct = sum(case[m] for m in members)
        if case_ct and case_ct != len(members):
            sibship_groups.append(members)

    unrelated = [iid for iid in kept
                 if iid not in in_family and iid not in in_sibship
                 and founder[iid]]
    if not (any(case[u] for u in unrelated) and
            any(not case[u] for u in unrelated)):
        unrelated = []

    with open(vcf_path) as f:
        for line in f:
            if line.startswith('##'):
                continue
            vcf_samples = line.split()[9:]
            break
        print('#ID\tOBS_CT\tEXP_CT\tCHISQ')
        for line in f:
            g = line.rstrip('\n').split('\t')
            if g[0] in ('X', 'Y', 'MT', 'M', '23', '24', '26'):
                # --dfam has no haploid/chrMT implementation.
                continue
            vid = g[2]
            geno = {}
            for iid, field in zip(vcf_samples, g[9:]):
                call = field.split(':')[0]
                if '.' in call:
                    geno[iid] = 3
                else:
                    a, b = call.replace('|', '/').split('/')
                    geno[iid] = int(a) + int(b)
            if counted_allele.get(vid, g[4]) == g[3]:
                for iid in geno:
                    if geno[iid] != 3:
                        geno[iid] = 2 - geno[iid]
            # Mendel errors are read from the unmodified genotypes.
            orig = dict(geno)
            for (pat, mat), children in families.items():
                for child in children:
                    if orig[child] == 3:
                        continue
                    err = MENDEL_TABLE[orig[pat] + orig[mat] * 4 +
                                       orig[child] * 16]
                    if not err:
                        continue
                    geno[child] = 3
                    if err & 0x100:
                        geno[pat] = 3
                    if err & 0x10000:
                        geno[mat] = 3
            twice_numer = 0
            quad_denom = 0
            total_ct = 0
            twice_total_expected = 0
            numer = 0.0
            denom = 0.0
            total_expected = 0.0
            for pat, mat, children in all_case_groups:
                p = PARENTAL.get((geno[pat], geno[mat]), 0)
                if not p:
                    continue
                case_ct = 0
                case_alt = 0
                for child in children:
                    if geno[child] == 3:
                        continue
                    case_ct += 1
                    case_alt += geno[child]
                if case_ct:
                    twice_numer += 2 * case_alt - case_ct * p
                    quad_denom += (2 - (p & 1)) * case_ct
                    total_ct += case_alt
                    twice_total_expected += case_ct * p
            for pat, mat, children in mixed_groups:
                p = PARENTAL.get((geno[pat], geno[mat]), 0)
                cc = ch = ce = uc = uh = ue = 0
                for child in children:
                    cg = geno[child]
                    if cg == 3:
                        continue
                    if case[child]:
                        cc += 1
                        ch += (cg == 2)
                        ce += (cg == 1)
                    else:
                        uc += 1
                        uh += (cg == 2)
                        ue += (cg == 1)
                if not cc:
                    continue
                if not p:
                    t, n, d, e = sibship_calc(cc, ch, ce, uc, uh, ue, plink19)
                    total_ct += t
                    numer += n
                    denom += d
                    total_expected += e
                else:
                    case_alt = 2 * ch + ce
                    twice_numer += 2 * case_alt - cc * p
                    quad_denom += (2 - (p & 1)) * (cc + uc)
                    total_ct += case_alt
                    twice_total_expected += cc * p
            numer += 0.5 * twice_numer
            denom += 0.25 * quad_denom
            total_expected += 0.5 * twice_total_expected
            for members in sibship_groups:
                cc = ch = ce = uc = uh = ue = 0
                for iid in members:
                    cg = geno[iid]
                    if cg == 3:
                        continue
                    if case[iid]:
                        cc += 1
                        ch += (cg == 2)
                        ce += (cg == 1)
                    else:
                        uc += 1
                        uh += (cg == 2)
                        ue += (cg == 1)
                if not cc:
                    continue
                t, n, d, e = sibship_calc(cc, ch, ce, uc, uh, ue, plink19)
                total_ct += t
                numer += n
                denom += d
                total_expected += e
            if unrelated:
                cc = ch = ce = uc = uh = ue = 0
                for iid in unrelated:
                    cg = geno[iid]
                    if cg == 3:
                        continue
                    if case[iid]:
                        cc += 1
                        ch += (cg == 2)
                        ce += (cg == 1)
                    else:
                        uc += 1
                        uh += (cg == 2)
                        ue += (cg == 1)
                case_alt = 2 * ch + ce
                hom = ch + uh
                het = ce + ue
                n_obs = cc + uc
                if not (n_obs <= 1 or ((not hom) and (not het)) or
                        hom == n_obs or het == n_obs):
                    total_ct += case_alt
                    if (not cc) or (not uc):
                        total_expected += case_alt
                    else:
                        p = cc / n_obs
                        cluster_alt = 2 * hom + het
                        expected = p * cluster_alt
                        var = (expected * (2 * n_obs - cluster_alt) * uc /
                               (n_obs * (2 * n_obs - 1)))
                        numer += case_alt - expected
                        denom += var
                        total_expected += expected
            chisq = 'NA' if denom == 0.0 else '%.10g' % (numer * numer / denom)
            print('%s\t%d\t%.10g\t%s' % (vid, total_ct, total_expected, chisq))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Independent least-squares oracle for --epistasis.

Reads a PLINK 1 fileset and an optional covariate file, fits
  phenotype ~ 1 + A + B + A*B + covariates
over the samples with both genotypes called, and prints the interaction term's
coefficient, standard error and t-statistic for every pair.
"""
import sys
import numpy as np

prefix = sys.argv[1]
covar_fname = sys.argv[2] if len(sys.argv) > 2 else None

fam = [l.split() for l in open(prefix + '.fam') if l.strip()]
sample_ids = [(f[0], f[1]) for f in fam]
pheno = np.array([float(f[5]) for f in fam])
n = len(fam)
bim = [l.split() for l in open(prefix + '.bim') if l.strip()]
ids = [b[1] for b in bim]

raw = open(prefix + '.bed', 'rb').read()
assert raw[:3] == b'\x6c\x1b\x01'
stride = (n + 3) // 4
codes = np.frombuffer(raw[3:], dtype=np.uint8).reshape(len(bim), stride)
bits = np.unpackbits(codes, axis=1, bitorder='little').reshape(len(bim), stride * 4, 2)
vals = bits[:, :n, 0].astype(np.int16) + 2 * bits[:, :n, 1].astype(np.int16)
# plink1 code 0 = hom A1, 1 = missing, 2 = het, 3 = hom A2; plink2 counts A1
geno = np.full(vals.shape, -1, dtype=np.float64)
geno[vals == 0] = 2.0
geno[vals == 2] = 1.0
geno[vals == 3] = 0.0
missing = (vals == 1)

keep = np.isfinite(pheno)
covars = np.zeros((n, 0))
if covar_fname:
    header = open(covar_fname).readline().split()
    rows = {}
    for line in open(covar_fname).readlines()[1:]:
        f = line.split()
        rows[(f[0], f[1])] = [float('nan') if x in ('NA', 'nan', '-9') else float(x) for x in f[2:]]
    covars = np.array([rows.get(sid, [np.nan] * (len(header) - 2)) for sid in sample_ids])
    keep &= np.isfinite(covars).all(axis=1)

print('ID1\tID2\tBETA_INT\tSE\tT_STAT')
out = []
for i in range(len(bim)):
    for j in range(i + 1, len(bim)):
        sel = keep & (~missing[i]) & (~missing[j])
        a = geno[i][sel]
        b = geno[j][sel]
        x = np.column_stack([np.ones(sel.sum()), a, b, a * b, covars[sel]])
        y = pheno[sel]
        xtx = x.T @ x
        try:
            inv = np.linalg.inv(xtx)
        except np.linalg.LinAlgError:
            continue
        beta = inv @ (x.T @ y)
        resid = y - x @ beta
        df = sel.sum() - x.shape[1]
        if df <= 0:
            continue
        sigma2 = (resid @ resid) / df
        se = np.sqrt(inv[3, 3] * sigma2)
        out.append('%s\t%s\t%.10g\t%.10g\t%.10g' % (ids[i], ids[j], beta[3], se, beta[3] / se))
print('\n'.join(out))

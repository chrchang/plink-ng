#!/usr/bin/env python3
"""Checks the INFO/AC and INFO/AN fields written by --export vcf vcf-info=.

Usage: check_ac_an.py <.vcf> <requested keys, e.g. AC,AN>

Every requested key must appear exactly once per variant (INFO/AC is omitted
when ALT is '.'), with values equal to the allele counts recomputed from the
GT field as written, i.e. with the ploidy of each exported call.  Keys that
were not requested must not have been touched, so they aren't checked.  The
header must declare each requested key exactly once, with the standard
Number/Type.  Standard library only.
"""

import sys

HEADER_LINES = {
    "AC": '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in '
    'genotypes, for each ALT allele, in the same order as listed">',
    "AN": '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of '
    'alleles in called genotypes">',
}


def fail(msg):
    sys.stderr.write("check_ac_an.py: " + msg + "\n")
    sys.exit(1)


def main():
    if len(sys.argv) != 3:
        fail("usage: check_ac_an.py <.vcf> <keys>")
    keys = sys.argv[2].split(",")
    header_ct = {key: 0 for key in HEADER_LINES}
    variant_ct = 0
    with open(sys.argv[1]) as vcf:
        for line in vcf:
            line = line.rstrip("\n")
            if line.startswith("##"):
                for key in HEADER_LINES:
                    if line.startswith("##INFO=<ID=" + key + ","):
                        header_ct[key] += 1
                        if key in keys and line != HEADER_LINES[key]:
                            fail("unexpected header line: " + line)
                continue
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            alts = [] if fields[4] == "." else fields[4].split(",")
            fmt = fields[8].split(":")
            if fmt[0] != "GT":
                fail("GT is not the first FORMAT field")
            allele_cts = [0] * (len(alts) + 1)
            an = 0
            for sample_field in fields[9:]:
                gt = sample_field.split(":")[0]
                for allele in gt.replace("|", "/").split("/"):
                    if allele != ".":
                        allele_cts[int(allele)] += 1
                        an += 1
            info = {}
            for entry in fields[7].split(";"):
                key, _, val = entry.partition("=")
                if key in info:
                    fail("duplicate INFO/%s for variant %s" % (key, fields[2]))
                info[key] = val
            expected = {"AN": str(an)}
            if alts:
                expected["AC"] = ",".join(str(ct) for ct in allele_cts[1:])
            for key in keys:
                if key not in expected:
                    if key in info:
                        fail("INFO/%s present for variant %s with no ALT allele"
                             % (key, fields[2]))
                    continue
                if info.get(key) != expected[key]:
                    fail("variant %s: INFO/%s=%s, expected %s" %
                         (fields[2], key, info.get(key), expected[key]))
            variant_ct += 1
    for key in keys:
        if header_ct[key] != 1:
            fail("%d ##INFO=<ID=%s header lines" % (header_ct[key], key))
    if not variant_ct:
        fail("no variants")


if __name__ == "__main__":
    main()

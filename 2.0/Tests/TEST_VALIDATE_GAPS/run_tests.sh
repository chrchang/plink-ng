#!/bin/bash

# --validate checks added after the review against pgen_spec:
# - dosage/hardcall consistency (hardcall missing, or within 0.5 of the total
#   ALT dosage), including "missing unconditional dosage => missing hardcall";
# - phased-dosage delta -32768 (which used to slip through uint16 wraparound
#   when the dosage was 0 or 2);
# - fixed-width dosage storage modes 0x03/0x04 (previously not validated);
# - exact layout: no gap before the first variant record, no bytes after the
#   last record/footer, .pgi ends after its header; header/footer extension
#   structure; padding (all-zero) and UTF-8 writer-ID header extensions;
# - multiallelic hardcall track with no mode byte.

set -exo pipefail

python3 make_pgens.py

for f in tmp_ok_dosage tmp_ok_ext tmp_ok_ext_pgi tmp_ok_pgi tmp_ok_fw3 tmp_ok_fw4; do
    $1/plink2 $2 $3 --pfile $f --validate --out ${f}_v
done

# Extensions and layout don't change the decoded data.
for f in tmp_ok_dosage tmp_ok_ext tmp_ok_ext_pgi tmp_ok_pgi; do
    $1/plink2 $2 $3 --pfile $f --export vcf vcf-dosage=HDS --out ${f}_e
done
for f in tmp_ok_ext tmp_ok_ext_pgi tmp_ok_pgi; do
    diff <(grep -v '^##' tmp_ok_dosage_e.vcf) <(grep -v '^##' ${f}_e.vcf)
done

expect_reject() {
    if $1/plink2 $2 $3 --pfile tmp_bad_$4 --validate --out tmp_bad_$4_v; then
        echo "expected --validate to reject tmp_bad_$4"
        exit 1
    fi
    # the log is word-wrapped
    tr '\n' ' ' < tmp_bad_$4_v.log | grep -qF "$5"
}

D="$1"
A="$2"
B="$3"
expect_reject "$D" "$A" "$B" incons_uncond "is more than 0.5 away from its (non-missing) hardcall"
expect_reject "$D" "$A" "$B" incons_uncond_hi "is more than 0.5 away from its (non-missing) hardcall"
expect_reject "$D" "$A" "$B" incons_list "is more than 0.5 away from its (non-missing) hardcall"
expect_reject "$D" "$A" "$B" incons_bits "is more than 0.5 away from its (non-missing) hardcall"
expect_reject "$D" "$A" "$B" missing_uncond "is missing, but its hardcall is not"
expect_reject "$D" "$A" "$B" range_uncond "dosage is greater than 2"
expect_reject "$D" "$A" "$B" dphase_m32768_hi "Invalid unconditional phased-dosages"
expect_reject "$D" "$A" "$B" dphase_m32768_lo "Invalid unconditional phased-dosages"
expect_reject "$D" "$A" "$B" dphase_m32768_list "Invalid phased-dosage track"
expect_reject "$D" "$A" "$B" mh_no_mode_byte "is missing its mode byte"
expect_reject "$D" "$A" "$B" gap_header "does not start immediately after the header ("
expect_reject "$D" "$A" "$B" gap_header_ext "does not start immediately after the header extensions"
expect_reject "$D" "$A" "$B" gap_pgi "does not start immediately after the 3-byte file header"
expect_reject "$D" "$A" "$B" pgi_trailing "pgi file size is inconsistent with its contents"
expect_reject "$D" "$A" "$B" trailing "extra bytes after the last variant record"
expect_reject "$D" "$A" "$B" trailing_footer "extra bytes after its footer"
expect_reject "$D" "$A" "$B" footer_gap "footer does not start immediately after the last variant record"
expect_reject "$D" "$A" "$B" padding_nonzero "padding header extension contains a nonzero byte"
expect_reject "$D" "$A" "$B" writer_id_utf8 "writer-identifier header extension is not valid UTF-8"
expect_reject "$D" "$A" "$B" writer_id_truncated "writer-identifier header extension is not valid UTF-8"
expect_reject "$D" "$A" "$B" flag_bit_256 "has a bit past 255 set"
expect_reject "$D" "$A" "$B" fw3_incons "is more than 0.5 away from its (non-missing) hardcall"
expect_reject "$D" "$A" "$B" fw3_missing "is missing, but its hardcall is not"
expect_reject "$D" "$A" "$B" fw3_range "dosage is greater than 2"
expect_reject "$D" "$A" "$B" fw4_dphase "Invalid unconditional phased-dosages"

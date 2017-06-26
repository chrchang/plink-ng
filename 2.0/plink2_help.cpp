// This file is part of PLINK 2.00, copyright (C) 2005-2017 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


// necessary to include this instead of plink2_common so g_cmdline_format_str[]
// is known to have external linkage
#include "plink2_help.h"

#ifdef __cplusplus
namespace plink2 {
#endif

const char g_cmdline_format_str[] = "\n  plink2 [input flag(s)...] {command flag(s)...} {other flag(s)...}\n  plink2 --help {flag name(s)...}\n\n";

uint32_t edit1_match(const char* s1, const char* s2, uint32_t len1, uint32_t len2) {
  // permit one difference of the following forms:
  // - inserted/deleted character
  // - replaced character
  // - adjacent pair of swapped characters
  uint32_t diff_found = 0;
  uint32_t pos = 0;
  if (len1 == len2) {
    while (pos < len1) {
      if (s1[pos] != s2[pos]) {
	if (diff_found) {
	  if ((diff_found == 2) || (s1[pos] != s2[pos - 1]) || (s1[pos - 1] != s2[pos])) {
	    return 0;
	  }
	}
	++diff_found;
      }
      ++pos;
    }
  } else if (len1 == len2 - 1) {
    do {
      if (s1[pos - diff_found] != s2[pos]) {
	if (diff_found) {
	  return 0;
	}
	++diff_found;
      }
      ++pos;
    } while (pos < len2);
  } else if (len1 == len2 + 1) {
    do {
      if (s1[pos] != s2[pos - diff_found]) {
	if (diff_found) {
	  return 0;
	}
	++diff_found;
      }
      ++pos;
    } while (pos < len1);
  } else {
    return 0;
  }
  return 1;
}

// for better or worse, when this is too small we find out quickly due to
// segfault...
CONSTU31(kMaxEqualHelpParams, 13);

typedef struct help_ctrl_struct {
  uint32_t iters_left;
  uint32_t param_ct;
  char** argv;
  uintptr_t unmatched_ct;
  uintptr_t* all_match_arr;
  uintptr_t* prefix_match_arr;
  uintptr_t* perfect_match_arr;
  uint32_t* param_slens;
  uint32_t preprint_newline;
} help_ctrl_t;

void help_print(const char* cur_params, help_ctrl_t* help_ctrl_ptr, uint32_t postprint_newline, const char* payload) {
  if (help_ctrl_ptr->param_ct) {
    strcpy(g_textbuf, cur_params);
    uint32_t cur_param_ct = 1;
    char* cur_param_start[kMaxEqualHelpParams];
    cur_param_start[0] = g_textbuf;
    char* textbuf_iter = strchr(g_textbuf, '\t');
    while (textbuf_iter) {
      *textbuf_iter++ = '\0';
      cur_param_start[cur_param_ct++] = textbuf_iter;
      textbuf_iter = strchr(textbuf_iter, '\t');
    }
    if (help_ctrl_ptr->iters_left) {
      const uint32_t orig_unmatched_ct = help_ctrl_ptr->unmatched_ct;
      if (help_ctrl_ptr->unmatched_ct) {
	uint32_t arg_uidx = 0;
	if (help_ctrl_ptr->iters_left == 2) {
	  for (uint32_t arg_idx = 0; arg_idx < orig_unmatched_ct; ++arg_idx, ++arg_uidx) {
	    arg_uidx = next_unset_unsafe(help_ctrl_ptr->all_match_arr, arg_uidx);
	    for (uint32_t cur_param_idx = 0; cur_param_idx < cur_param_ct; ++cur_param_idx) {
	      if (!strcmp(cur_param_start[cur_param_idx], help_ctrl_ptr->argv[arg_uidx])) {
		SET_BIT(arg_uidx, help_ctrl_ptr->perfect_match_arr);
		SET_BIT(arg_uidx, help_ctrl_ptr->prefix_match_arr);
		SET_BIT(arg_uidx, help_ctrl_ptr->all_match_arr);
		help_ctrl_ptr->unmatched_ct -= 1;
		break;
	      }
	    }
	  }
	} else {
          uint32_t cur_param_slens[kMaxEqualHelpParams];
	  for (uint32_t cur_param_idx = 0; cur_param_idx < cur_param_ct; ++cur_param_idx) {
	    cur_param_slens[cur_param_idx] = strlen(cur_param_start[cur_param_idx]);
	  }
	  for (uint32_t arg_idx = 0; arg_idx < orig_unmatched_ct; ++arg_idx, ++arg_uidx) {
	    arg_uidx = next_unset_unsafe(help_ctrl_ptr->all_match_arr, arg_uidx);
	    const uint32_t slen = help_ctrl_ptr->param_slens[arg_uidx];
	    for (uint32_t cur_param_idx = 0; cur_param_idx < cur_param_ct; ++cur_param_idx) {
	      if (cur_param_slens[cur_param_idx] > slen) {
		if (!memcmp(help_ctrl_ptr->argv[arg_uidx], cur_param_start[cur_param_idx], slen)) {
		  SET_BIT(arg_uidx, help_ctrl_ptr->prefix_match_arr);
		  SET_BIT(arg_uidx, help_ctrl_ptr->all_match_arr);
		  help_ctrl_ptr->unmatched_ct -= 1;
		  break;
		}
	      }
	    }
	  }
	}
      }
    } else {
      uint32_t cur_param_slens[kMaxEqualHelpParams];
      for (uint32_t cur_param_idx = 0; cur_param_idx < cur_param_ct; ++cur_param_idx) {
	cur_param_slens[cur_param_idx] = strlen(cur_param_start[cur_param_idx]);
      }
      uint32_t print_this = 0;
      for (uint32_t arg_uidx = 0; arg_uidx < help_ctrl_ptr->param_ct; ++arg_uidx) {
	if (IS_SET(help_ctrl_ptr->prefix_match_arr, arg_uidx)) {
	  if (!print_this) {
	    if (IS_SET(help_ctrl_ptr->perfect_match_arr, arg_uidx)) {
	      for (uint32_t cur_param_idx = 0; cur_param_idx < cur_param_ct; ++cur_param_idx) {
		if (!strcmp(cur_param_start[cur_param_idx], help_ctrl_ptr->argv[arg_uidx])) {
		  print_this = 1;
		  break;
		}
	      }
	    } else {
	      const uint32_t slen = help_ctrl_ptr->param_slens[arg_uidx];
	      for (uint32_t cur_param_idx = 0; cur_param_idx < cur_param_ct; ++cur_param_idx) {
		if (cur_param_slens[cur_param_idx] > slen) {
		  if (!memcmp(help_ctrl_ptr->argv[arg_uidx], cur_param_start[cur_param_idx], slen)) {
		    print_this = 1;
		    break;
		  }
		}
	      }
	    }
	  }
	} else {
	  for (uint32_t cur_param_idx = 0; cur_param_idx < cur_param_ct; ++cur_param_idx) {
	    if (edit1_match(cur_param_start[cur_param_idx], help_ctrl_ptr->argv[arg_uidx], cur_param_slens[cur_param_idx], help_ctrl_ptr->param_slens[arg_uidx])) {
	      print_this = 1;
	      if (!IS_SET(help_ctrl_ptr->all_match_arr, arg_uidx)) {
		SET_BIT(arg_uidx, help_ctrl_ptr->all_match_arr);
		help_ctrl_ptr->unmatched_ct -= 1;
	      }
	      break;
	    }
	  }
	}
      }
      if (print_this) {
	const uint32_t payload_slen = strlen(payload);
	const char* payload_end;
	if (payload[payload_slen - 2] == '\n') {
	  payload_end = &(payload[payload_slen - 1]);
	} else {
	  payload_end = &(payload[payload_slen]);
	}
	if (help_ctrl_ptr->preprint_newline) {
	  putc_unlocked('\n', stdout);
	}
	help_ctrl_ptr->preprint_newline = postprint_newline;
	const char* payload_iter = payload;
	do {
	  const char* line_end = (const char*)rawmemchr(payload_iter, '\n') + 1;
	  uint32_t line_slen = (uint32_t)(line_end - payload_iter);
	  if (line_slen > 2) {
	    payload_iter = &(payload_iter[2]);
	    line_slen -= 2;
	  }
	  memcpyx(g_textbuf, payload_iter, line_slen, 0);
	  fputs(g_textbuf, stdout);
	  payload_iter = line_end;
	} while (payload_iter < payload_end);
      }
    }
  } else {
    fputs(payload, stdout);
  }
}

pglerr_t disp_help(uint32_t param_ct, char** argv) {
  // yes, this is overkill.  But it should be a good template for other
  // command-line programs to use.
  uint32_t param_ctl = BITCT_TO_WORDCT(param_ct);
  pglerr_t reterr = kPglRetSuccess;
  help_ctrl_t help_ctrl;
  uint32_t arg_uidx;
  uint32_t arg_idx;
  uint32_t net_unmatched_ct;
  int32_t col_num;
  int32_t leading_dashes;
  help_ctrl.iters_left = param_ct? 2 : 0;
  help_ctrl.param_ct = param_ct;
  help_ctrl.argv = argv;
  help_ctrl.unmatched_ct = param_ct;
  help_ctrl.param_slens = nullptr;
  help_ctrl.all_match_arr = nullptr;
  help_ctrl.argv = nullptr;
  if (param_ct) {
    if (pgl_malloc(param_ct * sizeof(int32_t), &help_ctrl.param_slens) ||
	pgl_malloc(param_ctl * 3 * sizeof(intptr_t), &help_ctrl.all_match_arr)) {
      goto disp_help_ret_NOMEM;
    }
    leading_dashes = 0;
    for (arg_uidx = 0; arg_uidx < param_ct; arg_uidx++) {
      if (argv[arg_uidx][0] == '-') {
	leading_dashes = 1;
	break;
      }
    }
    if (leading_dashes) {
      if (pgl_malloc(param_ct * sizeof(intptr_t), &help_ctrl.argv)) {
	goto disp_help_ret_NOMEM;
      }
      for (arg_uidx = 0; arg_uidx < param_ct; arg_uidx++) {
	if (argv[arg_uidx][0] == '-') {
	  if (argv[arg_uidx][1] == '-') {
	    help_ctrl.argv[arg_uidx] = &(argv[arg_uidx][2]);
	  } else {
	    help_ctrl.argv[arg_uidx] = &(argv[arg_uidx][1]);
	  }
	} else {
	  help_ctrl.argv[arg_uidx] = argv[arg_uidx];
	}
      }
    } else {
      help_ctrl.argv = argv;
    }
    for (arg_idx = 0; arg_idx < param_ct; arg_idx++) {
      help_ctrl.param_slens[arg_idx] = strlen(help_ctrl.argv[arg_idx]);
    }
    fill_ulong_zero(param_ctl * 3, help_ctrl.all_match_arr);
    help_ctrl.prefix_match_arr = &(help_ctrl.all_match_arr[param_ctl]);
    help_ctrl.perfect_match_arr = &(help_ctrl.all_match_arr[param_ctl * 2]);
    help_ctrl.preprint_newline = 1;
  } else {
    help_ctrl.argv = nullptr;
    fputs(
"\nIn the command line flag definitions that follow,\n"
"  * [square brackets] denote a required parameter, where the text between the\n"
"    brackets describes its nature.\n"
"  * <angle brackets> denote an optional modifier (or if '|' is present, a set\n"
"    of mutually exclusive optional modifiers).  Use the EXACT text in the\n"
"    definition.\n"
"  * There's one exception to the angle brackets/exact text rule: when an angle\n"
"    bracket term ends with '=[value]', '[value]' designates a variable\n"
"    parameter.\n"
"  * {curly braces} denote an optional parameter, where the text between the\n"
"    braces describes its nature.\n"
"  * An ellipsis (...) indicates that you may enter multiple parameters of the\n"
"    specified type.\n"
"  * A \"column set descriptor\" is either\n"
"    1. a comma-separated sequence of column set names; this is interpreted as\n"
"       the full list of column sets to include.\n"
"    2. a comma-separated sequence of column set names, all preceded by '+' or\n"
"       '-'; this is interpreted as a list of changes to the default.\n"
, stdout);
    fputs(g_cmdline_format_str, stdout);
    fputs(
"Most " PROG_NAME_STR " runs require exactly one main input fileset.  The following flags\n"
"are available for defining its form and location:\n\n"
, stdout);
  }
  do {
    // explicit gzipped .pvar/.bim support was tried, and then rejected since
    // decompression was too slow
    // Zstd should have the necessary x86 performance characteristics, though
    help_print("pfile\tpgen\tbfile\tbed", &help_ctrl, 1,
"  --pfile [prefix] <vzs> : Specify .pgen + .pvar{.zst} + .psam prefix.\n"
"  --pgen [filename]      : Specify full name of .pgen/.bed file.\n"
	       );
    help_print("pfile\tpgen\tpvar\tpsam\tbfile\tbed\tbim\tfam\timport-dosage\tdosage", &help_ctrl, 1,
"  --pvar [filename]      : Specify full name of .pvar/.bim file.\n"
"  --psam [filename]      : Specify full name of .psam/.fam file.\n\n"
	       );
    help_print("bfile\tbpfile\tbed\tbim\tfam", &help_ctrl, 1,
"  --bfile [prefix] <vzs> : Specify .bed + .bim{.zst} + .fam prefix.\n"
"  --bpfile [prefx] <vzs> : Specify .pgen + .bim{.zst} + .fam prefix.\n\n"
	       );
    help_print("vcf\tbcf\tkeep-autoconv", &help_ctrl, 1,
"  --keep-autoconv    : When importing non-PLINK-binary data, don't delete\n"
"                       autogenerated binary fileset at end of run.\n\n"
	       );
    help_print("bfile\tfam", &help_ctrl, 1,
"  --no-fid           : .fam file does not contain column 1 (family ID).\n"
"  --no-parents       : .fam file does not contain columns 3-4 (parents).\n"
"  --no-sex           : .fam file does not contain column 5 (sex).\n"
"  --no-pheno         : .fam file does not contain column 6 (phenotype).\n\n"
	       );
    // todo: allele fraction import.  but first need to see how it's
    // represented in practice, since it isn't in the spec...
    help_print("vcf\tbcf\tpsam", &help_ctrl, 1,
"  --vcf [filename] <dosage=[field]>\n"
"  --bcf [filename] <dosage=[field]>  (not implemented yet)\n"
"    Specify full name of .vcf{.gz|.zst} or BCF2 file to import.\n"
"    * These can be used with --psam.\n"
"    * By default, dosage information is not imported.  To import the GP field\n"
"      (must be VCFv4.3-style 0..1, one probability per possible genotype), add\n"
"      'dosage=GP'.  'dosage=DS' (or anything else) causes the named field to be\n"
"      interpreted as a Minimac3-style dosage.\n\n"
	       );
    help_print("data\tgen\tbgen\tsample\thaps\tlegend", &help_ctrl, 1,
"  --data [filename prefix] <ref-first | ref-second> <gzs>\n"
"  --bgen [filename] <snpid-chr> <ref-first | ref-second>\n"
"  --gen [filename] <ref-first | ref-second>\n"
"  --sample [filename]\n"
"    Specify an Oxford-format dataset to import.  --data specifies a .gen{.zst}\n"
"    + .sample pair, while --bgen specifies a BGEN v1.1+ file.\n"
"    * If a BGEN v1.2+ file contains sample IDs, it may be imported without a\n"
"      companion .sample file.\n"
"    * With 'snpid-chr', chromosome codes are read from the 'SNP ID' field\n"
"      instead of the usual chromosome field.\n"
"    * By default, the second allele for each variant is treated as a\n"
"      provisional reference allele.  To specify that the first (resp. second)\n"
"      allele really is always reference, add the 'ref-first' (resp.\n"
"      'ref-second') modifier.\n\n"
	       );
    // todo: make 'per' prefix modifiable
    help_print("haps\tlegend", &help_ctrl, 1,
"  --haps [filename] <ref-first | ref-second>\n"
"  --legend [filename] [chr code]\n"
"    Specify .haps {+ .legend} file(s) to import.\n"
"    * When --legend is specified, it's assumed that the --haps file doesn't\n"
"      contain header columns.\n"
"    * On chrX, the second male column may contain dummy '-' entries.  (However,\n"
"      PLINK currently cannot handle omitted male columns.)\n"
"    * If not used with --sample, new sample IDs are of the form 'per#/per#'.\n\n"
	       );
    help_print("map\timport-dosage\tdosage", &help_ctrl, 1,
"  --map [fname]      : Specify full name of .map file.\n"
	       );
    help_print("import-dosage\tdosage", &help_ctrl, 1,
"  --import-dosage [allele dosage file] <noheader> <skip0=[i]> <skip1=[j]>\n"
"                  <skip2=[k]> <dose1> <format=[m]> <ref-first | ref-second>\n"
"                  <single-chr=[code]> <chr-col-num=[#]> <pos-col-num=[#]>\n"
"    Specify PLINK 1.x-style dosage file to import.\n"
"    * You must also specify a companion .psam/.fam file.\n"
"    * By default, PLINK assumes that the file contains a header line, which has\n"
"      'SNP' in (1-based) column i+1, 'A1' in column i+j+2, 'A2' in column\n"
"      i+j+3, and sample FID/IIDs starting from column i+j+k+4.  (i/j/k are\n"
"      normally zero, but can be changed with 'skip0', 'skip1', and 'skip2'\n"
"      respectively.)  If such a header line is not present, use the 'noheader'\n"
"      modifier; samples will then be assumed to appear in the same order as\n"
"      they do in the .psam/.fam file.\n"
"    * You may specify a companion .map file.  If you do not,\n"
"      * 'single-chr=' can be used to specify that all variants are on the named\n"
"        chromosome.  Otherwise, you can use 'chr-col-num=' to read chromosome\n"
"        codes from the given (1-based) column number.\n"
"      * 'pos-col-num=' causes bp coordinates to be read from the given column\n"
"        number.\n"
"    * The 'format' modifier lets you specify the number of values used to\n"
"      represent each dosage.  'format=1' normally indicates a single 0..2 A1\n"
"      expected count; 'dose1' modifies this to a 0..1 frequency.  'format=2'\n"
"      (the default) indicates a 0..1 homozygous A1 likelihood followed by a\n"
"      0..1 het likelihood, while 'format=3' indicates 0..1 hom A1, 0..1 het,\n"
"      0..1 hom A2.\n\n"
	       );
    // todo: triallelic rate
    help_print("dummy", &help_ctrl, 1,
"  --dummy [sample ct] [SNP ct] {missing dosage freq} {missing pheno freq}\n"
"          <acgt | 1234 | 12> <pheno-ct=[count]> <scalar-pheno>\n"
"          <dosage-freq=[rate]>\n"
"    This generates a fake input dataset with the specified number of samples\n"
"    and SNPs.\n"
"    * By default, the missing dosage and phenotype frequencies are zero.\n"
"      These can be changed by providing 3rd and 4th numeric parameters.\n"
"    * By default, allele codes are As and Bs; this can be changed with the\n"
"      'acgt', '1234', or '12' modifier.\n"
"    * By default, one binary phenotype is generated.  'pheno-ct=' can be used\n"
"      to change the number of phenotypes, and 'scalar-pheno' causes these\n"
"      phenotypes to be normally distributed scalars.\n"
"    * By default, all (nonmissing) dosages are in {0,1,2}.  To make some of\n"
"      them take on decimal values, use 'dosage-freq='.  (These dosages are\n"
"      affected by --hard-call-threshold and --dosage-erase-threshold.)\n\n"
	       );
    if (!param_ct) {
      fputs(
"Output files have names of the form '" PROG_NAME_STR ".{extension}' by default.  You can\n"
"change the '" PROG_NAME_STR "' prefix with\n\n"
, stdout);
    }
    help_print("out", &help_ctrl, 1,
"  --out [prefix]     : Specify prefix for output files.\n\n"
	       );
    if (!param_ct) {
      fputs(
"Most runs also require at least one of the following commands:\n\n"
, stdout);
    }
    help_print("make-pgen\tmake-bpgen\tmake-bed\tmake-just-pvar\tmake-just-psam", &help_ctrl, 1,
"  --make-pgen <vzs> <format=[code]> <trim-alts>\n"
"              <erase-phase> <erase-dosage>\n"
"              <pvar-cols=[col set descriptor]> <psam-cols=[col set descriptor]>\n"
"  --make-bpgen <vzs> <format=[code]> <trim-alts>\n"
"               <erase-phase> <erase-dosage>\n"
"  --make-bed <vzs> <trim-alts>\n"
	       /*
"  --make-pgen <vzs> <format=[code]> <multiallelics=[mode]> <trim-alts>\n"
"              <erase-alt2+> <erase-phase> <erase-dosage>\n"
"              <pvar-cols=[col set descriptor]> <psam-cols=[col set descriptor]>\n"
"  --make-bpgen <vzs> <format=[code]> <multiallelics=[mode]> <trim-alts>\n"
"               <erase-alt2+> <erase-phase> <erase-dosage>\n"
"  --make-bed <vzs> <multiallelics=[split mode]> <trim-alts>\n"
	       */
"    Create a new PLINK binary fileset (--make-pgen = .pgen + .pvar{.zst} +\n"
"    .psam, --make-bpgen = .pgen + .bim{.zst} + .fam).\n"
"    * Unlike the automatic text-to-binary converters (which only heed\n"
"      chromosome filters), this supports all of " PROG_NAME_STR "'s filtering flags.\n"
"    * The 'vzs' modifier causes the variant file (.pvar/.bim) to be\n"
"      Zstd-compressed.\n"
"    * The 'format' modifier requests an uncompressed fixed-variant-width .pgen\n"
"      file.  (These do not directly support multiallelic variants.)  The\n"
"      following format code is currently supported:\n"
"        2: just like .bed, except with an extended (12-byte instead of 3-byte)\n"
"           header containing variant/sample counts, and rotated genotype codes\n"
"           (00 = hom ref, 01 = het, 10 = hom alt, 11 = missing).\n"
	       /*
"        3: unphased dosage data\n"
"        4: phased dosage data\n"
	       */
    // Commented out since, while this is on the roadmap, it isn't implemented
    // yet.  (This also applies to other commented-out help text.)
	       /*
"    * The 'multiallelics' modifier (alias: 'm') specifies a merge or split\n"
"      mode.  The following modes are currently supported (well, not yet):\n"
"      * '-': Split all multiallelic records.\n"
"      * '-snps': Split SNP-only multiallelic records.\n"
"      * '+'/'+both': Adjacent variants with identical CHROM/POS/REF are\n"
"                     classified as SNPs and non-SNPs; SNPs are merged into one\n"
"                     variant, and non-SNPs are merged into another.\n"
"      * '+snps': Similar to '+both', except only SNPs are merged.\n"
"      * '+any': All adjacent biallelic variants with identical CHROM/POS/REF\n"
"                are merged into a single multiallelic variant.\n"
"      If a variant ID template was specified with --set-[missing/all]-var-ids,\n"
"      it is applied to the newly created variants.  Otherwise, the ID is set to\n"
"      the missing value.\n"
"      When merging, the new variant gets the lowest QUAL and the union of the\n"
"      FILTER values.\n"
"      INFO splitting/merging and left-alignment and normalization of indels are\n"
"      not currently supported.  'bcftools norm' (possibly on a single-sample\n"
"      file) can be used for this.\n"
"    * The 'trim-alts' modifier causes alternate alleles not present in the\n"
"      dataset after filtering to be removed.\n"
"    * The 'erase-alt2+' modifier causes alt alleles past the first to be\n"
"      removed; affected genotypes are set to missing.  (trim-alts happens\n"
"      first.)\n"
	       */
"    * The 'erase-phase' and 'erase-dosage' modifiers prevent phase and dosage\n"
"      information from being written to the new .pgen.\n"
	       /*
"    * When the 'multiallelics=', 'trim-alts', and/or 'erase-...' modifier is\n"
"      present, --make-bed/--make-{b}pgen cannot be combined with other\n"
"      commands.  (They can be combined with other filters.)\n"
	       */
"    * The first five columns of a .pvar file are always #CHROM/POS/ID/REF/ALT.\n"
"      Supported optional .pvar column sets are:\n"
"        xheader: All ## header lines (yeah, this is technically not a column).\n"
"                 Without this, only the #CHROM header line is kept.\n"
"        maybequal: QUAL.  Omitted if all loaded values are missing.\n"
"        qual: Force QUAL column to be written even when empty.\n"
"        maybefilter: FILTER.  Omitted if all loaded values are missing.\n"
"        filter: Force FILTER column to be written even when empty.\n"
"        maybeinfo: INFO.  Omitted if all loaded values are missing, or if\n"
"                   INFO:PR is the only subfield.\n"
"        info: Force INFO column to be written.\n"
"        maybecm: Centimorgan coordinate.  Omitted if all loaded values are 0.\n"
"        cm: Force CM column to be written even when empty.\n"
"      The default is xheader,maybequal,maybefilter,maybeinfo,maybecm.\n"
"    * The first two columns of a .psam file are always #FID/IID.  Supported\n"
"      optional .psam column sets are:\n"
"        maybesid: Sample disambiguation ID (useful when multiple samples are\n"
"                  collected from a single organism), '0' = missing.  Omitted if\n"
"                  all loaded values are missing.\n"
"        sid: Force SID column to be written even when empty.\n"
"        maybeparents: Father and mother IIDs, '0' = missing.  Omitted if all\n"
"                      loaded values are missing.\n"
"        parents: Force PAT and MAT columns to be written even when empty.\n"
"        sex: '1'/'M'/'m' = male, '2'/'F'/'f' = female, 'NA'/'0' = missing.\n"
"        pheno1: First active phenotype.  If none, all column entries are set to\n"
"                the --output-missing-phenotype string.\n"
"        phenos: All active phenotypes, if any.  (Can be combined with pheno1 to\n"
"                force at least one phenotype column to be written.)\n"
"      The default is maybesid,maybeparents,sex,phenos.\n\n"
	       );
    help_print("make-just-pvar\tmake-just-psam\tmake-just-bim\tmake-just-fam\twrite-cluster\n", &help_ctrl, 1,
"  --make-just-pvar <zs> <cols=[column set descriptor]>\n"
"  --make-just-psam <cols=[column set descriptor]>\n"
"  --make-just-bim <zs>\n"
"  --make-just-fam\n"
"    Variants of --make-pgen/--make-bed which only write a new variant or sample\n"
"    file.  These don't always require an input genotype file.\n"
"    USE THESE CAUTIOUSLY.  It is very easy to desynchronize your binary\n"
"    genotype data and your sample/variant indexes if you use these commands\n"
"    improperly.  If you have any doubt, stick with --make-{b}pgen/--make-bed.\n\n"
	       );
    help_print("export\trecode", &help_ctrl, 1,
"  --export [output format(s)...] <01 | 12> <bgz> <id-delim=[char]>\n"
"    <id-paste=[column set descriptor]> <include-alt> <omit-nonmale-y> <spaces>\n"
"    <vcf-dosage=[field]> <ref-first> <bits=[#]>\n"
"    Create a new fileset with all filters applied.  The following output\n"
"    formats are supported:\n"
"    (actually, only A-transpose, bgen-1.1, ind-major-bed, haps, hapslegend,\n"
"    oxford, and vcf are implemented for now)\n"
"    * '23': 23andMe 4-column format.  This can only be used on a single\n"
"            sample's data (--keep may be handy), and does not support\n"
"            multicharacter allele codes.\n"
"    * 'A': Sample-major additive (0/1/2) coding, suitable for loading from R.\n"
"           If you need uncounted alleles to be named in the header line, add\n"
"           the 'include-alt' modifier.\n"
"    * 'AD': Sample-major additive (0/1/2) + dominant (het=1/hom=0) coding.\n"
"            Also supports 'include-alt'.\n"
"    * 'A-transpose': Variant-major 0/1/2.\n"
"    * 'beagle': Unphased per-autosome .dat and .map files, readable by early\n"
"                BEAGLE versions.\n"
"    * 'beagle-nomap': Single .beagle.dat file.\n"
"    * 'bgen-1.x': Oxford-format .bgen + .sample.  For v1.2/v1.3, sample\n"
"                  identifiers are stored in the .bgen (with id-delim and\n"
"                  id-paste settings applied), and default precision is 16-bit\n"
"                  (use the 'bits' modifier to change this).\n"
"    * 'bimbam': Regular BIMBAM format.\n"
"    * 'bimbam-1chr': BIMBAM format, with a two-column .pos.txt file.  Does not\n"
"                     support multiple chromosomes.\n"
"    * 'fastphase': Per-chromosome fastPHASE files, with\n"
"                   .chr-[chr #].phase.inp filename extensions.\n"
"    * 'fastphase-1chr': Single .phase.inp file.  Does not support\n"
"                        multiple chromosomes.\n"
"    * 'haps', 'hapslegend': Oxford-format .haps + .sample{ + .legend}.  All\n"
"                            data must be biallelic and phased.\n"
"    * 'HV': Per-chromosome Haploview files, with .chr-[chr #][.ped + .info]\n"
"            filename extensions.\n"
"    * 'HV-1chr': Single Haploview .ped + .info file pair.  Does not support\n"
"                 multiple chromosomes.\n"
"    * 'ind-major-bed': PLINK 1 sample-major .bed (+ .bim + .fam).\n"
"    * 'lgen': PLINK 1 long-format (.lgen + .fam + .map), loadable with --lfile.\n"
"    * 'lgen-ref': .lgen + .fam + .map + .ref, loadable with --lfile +\n"
"                  --reference.\n"
"    * 'list': Single genotype-based list, up to 4 lines per variant.  To omit\n"
"              nonmale genotypes on the Y chromosome, add the 'omit-nonmale-y'\n"
"              modifier.\n"
"    * 'rlist': .rlist + .fam + .map fileset, where the .rlist file is a\n"
"                genotype-based list which omits the most common genotype for\n"
"                each variant.  Also supports 'omit-nonmale-y'.\n"
"    * 'oxford': Oxford-format .gen + .sample.\n"
"    * 'ped': PLINK 1 sample-major (.ped + .map), loadable with --file.\n"
"    * 'compound-genotypes': Same as 'ped', except that the space between each\n"
"                            pair of same-variant allele codes is removed.\n"
"    * 'structure': Structure-format.\n"
"    * 'transpose': PLINK 1 variant-major (.tped + .tfam), loadable with\n"
"                   --tfile.\n"
"    * 'vcf': VCFv4.3.  If PAR1 and PAR2 are present, they are automatically\n"
"             merged with chrX, with proper handling of chromosome codes and\n"
"             male ploidy.  If the 'bgz' modifier is added, the VCF file is\n"
"             block-gzipped.\n"
"             The 'id-paste' modifier controls which .psam columns are used to\n"
"             construct sample IDs (choices are fid, iid, maybesid, and sid;\n"
"             default is fid,iid,maybesid), while the 'id-delim' modifier sets\n"
"             the character between the ID pieces (default '_').\n"
"             By default, dosages are not exported; use 'vcf-dosage=GP' to\n"
"             export them as genotype posterior probabilities, or\n"
"             'vcf-dosage=DS' to export Minimac3-style dosages.\n"
	       // possible todo: pedigree output?
"    In addition,\n"
"    * When the output format only supports biallelic variants, multiallelic\n"
"      variants are downcoded to ref/alt1, not split.\n"
	       // todo: implement CPRA <-> CPR
"    * The '12' modifier causes alt1 alleles to be coded as '1' and ref alleles\n"
"      to be coded as '2', while '01' maps alt1 -> 0 and ref -> 1.\n"
"    * The 'spaces' modifier makes the output space-delimited instead of\n"
"      tab-delimited, whenever both are permitted.\n"
"    * For biallelic formats where it's unspecified whether the reference/major\n"
"      allele should appear first or second, --export defaults to second for\n"
"      compatibility with PLINK 1.9.  Use 'ref-first' to change this.\n\n"
	       );
    
    // don't bother with case/control or cluster-stratification any more, since
    // user can loop through subgroups and then use Unix cut/paste

    // todo: add optional columns for computed MAF (nothing here quite
    // corresponds to nonmajor_freqs when e.g. --maf-succ was specified) and
    // machr2 (this is probably the best home since, unlike --geno-counts and
    // --hardy, it's dosage-aware).
    help_print("freq\tmach-r2-filter", &help_ctrl, 1,
"  --freq <zs> <counts> <cols=[column set descriptor]> <bins-only>\n"
"         <refbins=[comma-separated bin boundaries] | refbins-file=[filename]>\n"
"         <alt1bins=[comma-separated bin boundaries] | alt1bins-file=[filename]>\n"
"    Empirical allele frequency report.  By default, only founders are\n"
"    considered.  Dosages are taken into account (e.g. heterozygous haploid\n"
"    calls count as 0.5).  chrM dosages are scaled to sum to 2.\n"
"    Supported column sets are:\n"
"      chrom: Chromosome ID.\n"
"      pos: Base-pair coordinate.\n"
"      (ID is always present, and positioned here.)\n"
"      ref: Reference allele.\n"
"      alt1: Alternate allele 1.\n"
"      alt: All alternate alleles, comma-separated.\n"
"      reffreq: Reference allele frequency/dosage.\n"
"      alt1freq: Alt1 frequency/dosage.\n"
"      altfreq: Comma-separated frequencies/dosages for all alternate alleles.\n"
"      freq: Similar to altfreq, except ref is also included at the start.\n"
"      eq: Comma-separated [allele]=[freq] for all present alleles.  (If no\n"
"          alleles are present, the column contains a single '.'.)\n"
"      eqz: Same as eq, except zero-counts are included.\n"
"      alteq/alteqz: Same as eq/eqz, except reference allele is omitted.\n"
"      numeq: 0=[freq],1=[freq], etc.  Zero-counts are omitted.\n"
"      altnumeq: Same as numeq, except reference allele is omitted.\n"
"      machr2: Empirical divided by theoretical variance quality metric.\n"
"      nobs: Number of allele observations.\n"
"    The default is chrom,ref,alt,altfreq,nobs.\n"
"    Additional .afreq.{ref,alt1}.bins (or .acount.{ref,alt1}.bins with\n"
"    'counts') file(s) are generated when 'refbins='/'refbins-file=' or\n"
"    'alt1bins='/'alt1bins-file=' is present; these report the total number of\n"
"    frequencies or counts in each left-closed, right-open interval.  (If you\n"
"    only want these histogram(s), and not the main report, add 'bins-only'.)\n\n"
	       );
    // this can't really handle dosages, so we specify "hardcall"
    help_print("geno-counts\tfreq\tfreqx\frqx", &help_ctrl, 1,
"  --geno-counts <zs> <cols=[column set descriptor]>\n"
"    Hardcall genotype count report (considering both alleles simultaneously in\n"
"    the diploid case).  Nonfounders are now included; use --keep-founders if\n"
"    this is a problem.  Heterozygous haploid calls are treated as missing.\n"
"    Supported column sets are:\n"
"      chrom: Chromosome ID.\n"
"      pos: Base-pair coordinate.\n"
"      (ID is always present, and positioned here.)\n"
"      ref: Reference allele.\n"
"      alt1: Alternate allele 1.\n"
"      alt: All alternate alleles, comma-separated.\n"
"      homref: Homozygous-ref count.\n"
"      refalt1: Heterozygous ref-alt1 count.\n"
"      refalt: Comma-separated het ref-altx counts.\n"
"      homalt1: Homozygous-alt1 count.\n"
"      altxy: Comma-separated altx-alty counts, in (1/1)-(1/2)-(2/2)-(1/3)-...\n"
"             order.\n"
"      xy: Similar to altxy, except the reference allele is treated as alt0,\n"
"          and the sequence starts (0/0)-(0/1)-(1/1)-(0/2)-...\n"
"      hapref: Haploid-ref count.\n"
"      hapalt1: Haploid-alt1 count.\n"
"      hapalt: Comma-separated haploid-altx counts.\n"
"      hap: Similar to hapalts, except ref is also included at the start.\n"
"      numeq: 0/0=[hom ref ct],0/1=[het ref-alt1],1/1=[hom alt1],...,0=[hap ref]\n"
"             etc.  Zero-counts are omitted.  (If all genotypes are missing, the\n"
"             column contains a single '.'.)\n"
"      missing: Number of missing genotypes.\n"
"      nobs: Number of (nonmissing) genotype observations.\n"
"    The default is chrom,ref,alt,homref,refalt,altxy,hapref,hapalt,missing.\n\n"
	       );
    // todo: add cluster-stratification
    help_print("missing", &help_ctrl, 1,
"  --missing <zs> <sample-only | variant-only> <scols=[column set descriptor]>\n"
"            <vcols=[column set descriptor]>\n"
"    Generate sample- and variant-based missing data reports (or just one report\n"
"    if 'sample-only'/'variant-only' is specified).\n"
"    Supported column sets in the sample-based report are:\n"
"      (FID and IID are always present, and positioned here.)\n"
"      maybesid: SID, if at least one nonmissing value is present.\n"
"      sid: Force SID column to be written even when empty.\n"
"      misspheno1: First active phenotype missing (Y/N)?  Always 'Y' if no\n"
"                  phenotypes are loaded.\n"
"      missphenos: A Y/N column for each loaded phenotype.  (Can be combined\n"
"                  with misspheno1 to force at least one such column.)\n"
"      nmissdosage: Number of missing dosages.\n"
"      nmiss: Number of missing hardcalls, not counting het haploids.\n"
"      nmisshh: Number of missing hardcalls, counting het haploids.\n"
"      hethap: Number of heterozygous haploid hardcalls.\n"
"      nobs: Denominator (male count on chrY, otherwise total sample count).\n"
"      fmissdosage: Missing dosage rate.\n"
"      fmiss: Missing hardcall rate, not counting het haploids.\n"
"      fmisshh: Missing hardcall rate, counting het haploids.\n"
"    The default is maybesid,missphenos,nmiss,nobs,fmiss.\n"
"    Supported column sets in the variant-based report are:\n"
"      chrom: Chromosome ID.\n"
"      pos: Base-pair coordinate.\n"
"      (ID is always present, and positioned here.)\n"
"      ref: Reference allele.\n"
"      alt1: Alternate allele 1.\n"
"      alt: All alternate alleles, comma-separated.\n"
"      nmissdosage: Number of missing dosages.\n"
"      nmiss: Number of missing hardcalls, not counting het haploids.\n"
"      nmisshh: Number of missing hardcalls, counting het haploids.\n"
"      hethap: Number of heterozygous haploid calls.\n"
"      nobs: Number of potentially valid calls.\n"
"      fmissdosage: Missing dosage rate.\n"
"      fmiss: Missing hardcall rate, not counting het haploids.\n"
"      fmisshh: Missing hardcall rate, counting het haploids.\n"
"      fhethap: Heterozygous haploid rate.\n"
"    The default is chrom,nmiss,nobs,fmiss.\n\n"
	       );
    help_print("hardy", &help_ctrl, 1,
"  --hardy <zs> <midp> <cols=[column set descriptor]>\n"
"    Hardy-Weinberg exact test p-value report(s).\n"
"    * For multiallelic variants, the test is based on the reference allele.\n"
"    * By default, only founders are considered; change this with --nonfounders.\n"
"    * chrX is now omitted from the main {output prefix}.hardy report.  Instead,\n"
"      (if present) it gets its own {output prefix}.hardy.x report based on the\n"
"      method described in Graffelman J, Weir BS (2016) Hardy-Weinberg\n"
"      equilibrium and the X chromosome.\n"
"    * There is currently no special handling of case/control phenotypes.\n"
"    Supported column sets are:\n"
"      chrom: Chromosome ID.\n"
"      pos: Base-pair coordinate.\n"
"      (ID is always present, and positioned here.)\n"
"      ref: Reference allele.\n"
"      alt1: Alternate allele 1.\n"
"      alt: All alternate alleles, comma-separated.\n"
"      gcounts: Hom-ref count, total number of ref-altx heterozygous calls, and\n"
"               total number of nonmissing calls with no reference allele.  On\n"
"               chrX, these are followed by male ref and male alt counts.\n"
"      gcount1col: gcounts values in a single comma-separated column.\n"
"      hetfreq: Observed and expected heterozygote frequencies.\n"
"      sexaf: Female and male ref allele frequencies (chrX only).\n"
"      femalep: Female-only p/midp-value (chrX only).\n"
"      p: Hardy-Weinberg equilibrium exact test p/midp-value.\n"
"    The default is chrom,ref,alt,gcounts,hetfreq,sexaf,p.\n\n"
	       );
    help_print("indep\tindep-pairwise", &help_ctrl, 1,
"  --indep-pairwise [window size]<kb> {step size (variant ct)} [r^2 threshold]\n"
"    Generate a list of variants in approximate linkage equilibrium.  With the\n"
"    'kb' modifier, the window size is in kilobase instead of variant count\n"
"    units.  (Pre-'kb' space is optional, i.e. '--indep-pairwise 500 kb 0.5' and\n"
"    and '--indep-pairwise 500kb 0.5' have the same effect.)\n"
"    The step size now defaults to 1 if it's unspecified, and *must* be 1 if the\n"
"    window is in kilobase units.\n"
"    Note that you need to rerun " PROG_NAME_STR " using --extract or --exclude on the\n"
"    .prune.in/.prune.out file to apply the list to another computation.\n\n"
	       );
    // todo: replace --indep-pairphase with method which takes fully phased
    // haplotypes as input (unphased het treated as missing).
    
    // for kinship estimation, LD pruning isn't really advisable (if more speed
    // is needed, the humble --bp-space may lead to a better approximation; and
    // in practice speed isn't an issue any more with --make-king)
    help_print("make-king\tmake-king-table", &help_ctrl, 1,
"  --make-king <square | square0 | triangle> <zs | bin | bin4>\n"
"    KING-robust kinship estimator, described by Manichaikul A, Mychaleckyj JC,\n"
"    Rich SS, Daly K, Sale M, Chen WM (2010) Robust relationship inference in\n"
"    genome-wide association studies.  By default, this writes a\n"
"    lower-triangular tab-delimited table of kinship coefficients to\n"
"    {output prefix}.king, and a list of the corresponding sample IDs to\n"
"    {output prefix}.king.id.  The first row of the .king file contains a single\n"
"    [genome 1-genome 2] kinship coefficient, the second row has the\n"
"    [genome 1-genome 3] and [genome 2-genome 3] kinship values in that order,\n"
"    etc.\n"
"    * Only autosomes are currently considered.\n"
"    * Pedigree information is currently ignored; the between-family estimator\n"
"      is used for all pairs.\n"
"    * If the 'square' or 'square0' modifier is present, a square matrix is\n"
"      written instead; 'square0' fills the upper right triangle with zeroes.\n"
"    * If the 'zs' modifier is present, the .king file is Zstd-compressed.\n"
"    * If the 'bin' modifier is present, a binary (square) matrix of\n"
"      double-precision floating point values, suitable for loading from R, is\n"
"      instead written to {output prefix}.king.bin.  ('bin4' specifies\n"
"      single-precision numbers instead.)  This can be combined with 'square0'\n"
"      if you still want the upper right zeroed out, or 'triangle' if you don't\n"
"      want to pad the upper right at all.\n"
"    * The computation can be subdivided with --parallel.\n"
"  --make-king-table <zs> <counts> <cols=[column set descriptor]>\n"
"    Similar to --make-king, except results are reported in the original .kin0\n"
"    text table format (with minor changes, e.g. row order is more friendly to\n"
"    incremental addition of samples), and --king-table-filter can be used to\n"
"    restrict the report to high kinship values.\n"
"    Supported column sets are:\n"
"      (FID and IID are always present, and positioned here.)\n"
"      maybesid: SID, if at least one nonmissing value is present.\n"
"      sid: Force SID column to be written even when empty.\n"
"      misspheno1: First active phenotype missing (Y/N)?  Always 'Y' if no\n"
"                  phenotypes are loaded.\n"
"      missphenos: A Y/N column for each loaded phenotype.  (Can be combined\n"
"                  with misspheno1 to force at least one such column.)\n"
"      id: FID1/ID1/FID2/ID2.\n"
"      maybesid: SID1/SID2, if at least one value is nonmissing.  Must be used\n"
"                with 'id'.\n"
"      sid: Force SID1/SID2 even when all values are missing.\n"
"      nsnp: Number of variants considered (autosomal, neither call missing).\n"
"      hethet: Proportion/count of considered call pairs which are het-het.\n"
"      ibs0: Proportion/count of considered call pairs which are opposite homs.\n"
"      ibs1: HET1_HOM2 and HET2_HOM1 proportions/counts.\n"
"      kinship: KING-robust between-family kinship estimator.\n"
"    The default is id,maybesid,nsnp,hethet,ibs0,kinship.  hethet/ibs0/ibs1\n"
"    values are proportions unless the 'counts' modifier is present.  If id is\n"
"    omitted, a .kin0.id file is also written.\n\n"
	       );
    help_print("make-rel\tmake-grm\tmake-grm-bin\tmake-grm-gz", &help_ctrl, 1,
"  --make-rel <cov> <meanimpute> <square | square0 | triangle> <zs | bin | bin4>\n"
"    Write a lower-triangular variance-standardized relationship matrix to\n"
"    {output prefix}.rel, and corresponding IDs to {output prefix}.rel.id.\n"
"    * It is usually best to perform this calculation on a variant set in\n"
"      approximate linkage equilibrium, with no very-low-MAF variants.\n"
"    * The 'cov' modifier removes the variance standardization step, causing a\n"
"      covariance matrix to be calculated instead.\n"
"    * The computation can be subdivided with --parallel.\n"
"  --make-grm-gz <cov> <meanimpute> <no-gz | zs>\n"
"  --make-grm-bin <cov> <meanimpute>\n"
"    --make-grm-gz causes the relationships to be written to GCTA's original\n"
"    gzipped list format, which describes one pair per line, while\n"
"    --make-grm-bin writes them in GCTA 1.1+'s single-precision triangular\n"
"    binary format.  Note that these formats explicitly report the number of\n"
"    valid observations (where neither sample has a missing call) for each pair,\n"
"    which is useful input for some scripts.\n\n"
	       );
#ifndef NOLAPACK
    // GRM, PCA, etc. based on major vs. nonmajor alleles
    // possible todo: have an 'approx2' mode which implements the flashpca 2.0
    //   algorithm, which does not require memory quadratic in the # of PCs
    help_print("pca", &help_ctrl, 1,
"  --pca {count} <approx | meanimpute> <sid>\n"
"  --pca var-wts {count} <approx | meanimpute> <sid> <vzs>\n"
"                <vcols=[col set descriptor]>\n"
"    Extracts top principal components from the variance-standardized\n"
"    relationship matrix.\n"
"    * It is usually best to perform this calculation on a variant set in\n"
"      approximate linkage equilibrium, with no very-low-MAF variants.\n"
"    * By default, 10 PCs are extracted; you can adjust this by passing a\n"
"      numeric parameter.  (Note that 10 is lower than the PLINK 1.9 default of\n"
"      20; this is due to the randomized algorithm's memory footprint growing\n"
"      quadratically w.r.t. the PC count.)\n"
"    * The 'approx' modifier causes the standard deterministic computation to be\n"
"      replaced with the randomized algorithm originally implemented for\n"
"      Galinsky KJ, Bhatia G, Loh PR, Georgiev S, Mukherjee S, Patterson NJ,\n"
"      Price AL (2016) Fast Principal-Component Analysis Reveals Convergent\n"
"      Evolution of ADH1B in Europe and East Asia.  This can be a good idea when\n"
"      you have >5k samples.\n"
"    * The randomized algorithm always uses mean imputation for missing genotype\n"
"      calls.  For comparison purposes, you can use the 'meanimpute' modifier to\n"
"      request this behavior for the standard computation.\n"
"    * The 'var-wts' modifier requests an additional .eigenvec.var file with PCs\n"
"      expressed as variant weights instead of sample weights.  When it's\n"
"      present, 'vzs' causes the .eigenvec.var file to be Zstd-compressed.\n"
"      'vcols' can be used to customize the report columns; supported column\n"
"      sets are:\n"
"        chrom: Chromosome ID.\n"
"        pos: Base-pair coordinate.\n"
"        (ID is always present, and positioned here.)\n"
"        ref: Reference allele.\n"
"        alt1: Alternate allele 1.\n"
"        alt: All alternate alleles, comma-separated.\n"
"        maj: Major allele.\n"
"        nonmaj: All nonmajor alleles, comma-separated.\n"
"        (PCs are always present, and positioned here.  Signs are w.r.t. the\n"
"        major, not necessarily reference, allele.)\n"
"      Default is chrom,maj,nonmaj.\n\n"
	       );
#endif
    help_print("king-cutoff\tmake-king\tmake-king-table\trel-cutoff\tgrm-cutoff", &help_ctrl, 1,
"  --king-cutoff {.king.bin + .king.id fileset prefix} [threshold]\n"
"    Exclude one member of each pair of samples with KING-robust kinship greater\n"
"    than the given threshold.  Remaining/excluded sample IDs are written to\n"
"    {output prefix}.king.cutoff.in + .king.cutoff.out.\n"
"    If present, the .king.bin file must be triangular (either precision is ok).\n\n"
	       );
    help_print("write-covar\twith-phenotype", &help_ctrl, 1,
"  --write-covar <cols=[column set descriptor]>\n"
"    If covariates are defined, an updated version (with all filters applied) is\n"
"    automatically written to {output prefix}.cov whenever --make-pgen,\n"
"    --make-just-psam, --export, or a similar command is present.  However, if\n"
"    you do not wish to simultaneously generate a new sample file, you can use\n"
"    --write-covar to just produce a pruned covariate file.\n"
"    Supported column sets are:\n"
"      maybesid: SID, if at least one nonmissing value is present.\n"
"      sid: Force SID column to be written even when empty.\n"
"      maybeparents: Father and mother IIDs, '0' = missing.  Omitted if all\n"
"                    loaded values are missing.\n"
"      parents: Force PAT and MAT columns to be written even when empty.\n"
"      sex: '1'/'M'/'m' = male, '2'/'F'/'f' = female, 'NA'/'0' = missing.\n"
"      pheno1: First active phenotype.  If none, all column entries are set to\n"
"              the --output-missing-phenotype string.\n"
"      phenos: All active phenotypes, if any.  (Can be combined with pheno1 to\n"
"              force at least one phenotype column to be written.)\n"
"      (Covariates are always present, and positioned here.)\n"
"    The default is just maybesid.\n\n"
	       );
    help_print("write-snplist", &help_ctrl, 1,
"  --write-snplist <zs>\n"
"    List all variants which pass your filters/inclusion thresholds.\n\n"
	       );
    help_print("glm\tlinear\tlogistic\tassoc", &help_ctrl, 1,
"  --glm <zs> <sex | no-x-sex> <genotypic | hethom | dominant | recessive>\n"
"        <interaction> <hide-covar> <intercept> <firth-fallback | firth>\n"
"        <cols=[col set descriptor]> <local-covar=[f]> <local-pvar=[f]>\n"
"        <local-psam=[f]> <local-omit-last | local-cats=[category ct]>\n"
	       // "        <perm | mperm=[value]> <perm-count>\n"
"    Basic association analysis on quantitative and/or case/control phenotypes.\n"
"    For each variant, a linear (for quantitative traits) or logistic (for\n"
"    case/control) regression is run with the phenotype as the dependent\n"
"    variable, and alt dosage and a constant-1 column as predictors.\n"
"    * For multiallelic variants, the total alt1 + alt2 + ... dosage is used.\n"
"    * By default, sex (male = 1, female = 2; note that this is a change from\n"
"      PLINK 1.x) is automatically added as a predictor for X chromosome\n"
"      variants, and no others.  The 'sex' modifier causes it to be added\n"
"      everywhere (except chrY), while 'no-x-sex' excludes it entirely.\n"
"    * The 'genotypic' modifier adds an additive effect/dominance deviation 2df\n"
"      joint test (0-2 and 0..1..0 coding), while 'hethom' uses 0..0..1 and\n"
"      0..1..0 coding instead.\n"
	       /*
"  If permutation is also requested, these\n"
"      modifiers cause permutation to be based on the joint test.\n"
	       */
"    * 'dominant' and 'recessive' specify a model assuming full dominance or\n"
"      recessiveness, respectively, for the ref allele.  I.e. the genotype\n"
"      column is recoded as 0..1..1 or 0..0..1, respectively.\n"
"    * 'interaction' adds genotype x covariate interactions to the model.\n"
	       /*
"  This\n"
"      cannot be combined with the usual permutation tests; use --tests to\n"
"      define the permutation test statistic instead.\n"
	       */
"    * Additional predictors can be added with --covar.  By default, association\n"
"      statistics are reported for all nonconstant predictors; 'hide-covar'\n"
"      suppresses covariate-only results, while 'intercept' causes intercepts\n"
"      to be reported.\n"
"    * For logistic regression, when the phenotype {quasi-}separates the\n"
"      genotype, an NA result will normally be reported.  To fall back on Firth\n"
"      logistic regression instead when the basic logistic regression fails to\n"
"      converge, add the 'firth-fallback' modifier.  To eliminate the special\n"
"      case and use Firth logistic regression everywhere, add 'firth'.\n"
"    * To add covariates which are not constant across all variants, add the\n"
"      'local-covar=', 'local-pvar=', and 'local-psam=' modifiers, and use full\n"
"      filenames for each.\n"
"      Normally, the local-covar file should have c * n real-valued columns,\n""      where the first c columns correspond to the first sample in the\n"
"      local-psam file, columns (c+1) to 2c correspond to the second sample,\n"
"      etc.; and the mth line correspond to the mth nonheader line of the\n"
"      local-pvar file.  (Variants outside of the local-pvar file are excluded\n"
"      from the regression.)  The local covariates are assigned the names\n"
"      LOCAL1, LOCAL2, etc.; to exclude the last local covariate from the\n"
"      regression (necessary if they are e.g. local ancestry coefficients which\n"
"      sum to 1), add 'local-omit-last'.\n"
"      Alternatively, with 'local-cats=[k]', the local-covar file is expected to\n"
"      have n columns with integer-valued entries in [1, k].  These category\n"
"      assignments are expanded into (k-1) local covariates in the usual manner.\n"
	       /*
"    * 'perm' normally causes an adaptive permutation test to be performed on\n"
"      the main effect, while 'mperm=[value]' starts a max(T) permutation test.\n"
"    * 'perm-count' causes the permutation test report to include counts instead\n"
"      of frequencies.\n"
	       */
// May want to change or leave out set-based test; punt for now.
"    The main report supports the following column sets:\n"
"      chrom: Chromosome ID.\n"
"      pos: Base-pair coordinate.\n"
"      (ID is always present, and positioned here.)\n"
"      ref: Reference allele.\n"
"      alt1: Alternate allele 1.\n"
"      alt: All alternate alleles, comma-separated.\n"
"      altcount: Alternate allele count (can be decimal with dosage data).\n"
"      totallele: Allele observation count (can be higher than --freq value, due\n"
"                 to inclusion of het haploids and chrX model).\n"
"      altcountcc: alt count in cases, then controls (case/control only).\n"
"      totallelecc: Case and control allele observation counts.\n"
"      altfreq: alt allele frequency.\n"
"      altfreqcc: alt frequency in cases, then controls (case/control only).\n"
"      machr2: Empirical divided by theoretical variance quality metric.\n"
"      firth: Reports whether Firth regression was used (firth-fallback only).\n"
"      test: Test identifier.  (Required unless only one test is run.)\n"
"      nobs: Number of samples in the regression.\n"
"      beta: Regression coefficient (for alternate allele).\n"
"      orbeta: Odds ratio for case/control, beta for quantitative traits.\n"
"      se: Standard error of beta/odds ratio.\n"
"      ci: Bounds of symmetric approximate confidence interval (requires --ci).\n"
"      t: T-statistic.\n"
"      p: Asymptotic p-value for t-statistic.\n"
"    The default is chrom,pos,ref,alt,firth,test,nobs,orbeta,se,ci,t,p.\n\n"
	       );
    help_print("score", &help_ctrl, 1,
"  --score [filename] {i} {j} {k} <header | header-read> <no-mean-imputation>\n"
"          <center | variance-standardize> <se> <zs>\n"
"          <list-variants | list-variants-zs> <cols=[col set descriptor]>\n"
"    Apply linear scoring system(s) to each sample.\n"
"    The input file should have one line per scored variant.  Variant IDs are\n"
"    read from column #i and allele codes are read from column #j, where i\n"
"    defaults to 1 and j defaults to i+1.\n"
"    * By default, a single column of input coefficients is read from column #k,\n"
"      where k defaults to j+1.  (--score-number can be used to specify multiple\n"
"      columns.)\n"
"    * The 'header' modifier causes the first nonempty line of the input file to\n"
"      be treated as an ignorable header line, while 'header-read' causes score\n"
"      column header(s) to be read and included in the report.\n"
"    * By default, copies of unnamed alleles contribute zero to score, while\n"
"      missing genotypes contribute an amount proportional to the loaded (via\n"
"      --read-freq) or imputed allele frequency.  To throw out missing\n"
"      observations instead (decreasing the denominator in the final average\n"
"      when this happens), use the 'no-mean-imputation' modifier.\n"
"    * You can use the 'center' modifier to shift all genotypes to mean zero, or\n"
"      'variance-standardize' to linearly transform the genotypes to mean-0,\n"
"      variance-1.  ('variance-standardize' cannot be used with chrX or MT.)\n"
"    * The 'se' modifier causes the score coefficients to be treated as\n"
"      independent standard errors; in this case, standard errors for the score\n"
"      average/sum are reported.  (Note that this will systematically\n"
"      underestimate standard errors when scored variants are in LD.)\n"
"    * The 'list-variants{-zs}' modifier causes variant IDs used for scoring to\n"
"      be written to [output prefix].sscore.vars{.zst}.\n"
"    The main report supports the following column sets:\n"
"      (FID and IID are always present, and positioned here.)\n"
"      maybesid: SID, if at least one nonmissing value is present.\n"
"      sid: Force SID column to be written even when empty.\n"
"      pheno1: First active phenotype.\n"
"      phenos: All active phenotypes, if any.\n"
"      nmissallele: Number of nonmissing alleles.\n"
"      denom: Denominator of score average (equal to nmissallele value when\n"
"             'no-mean-imputation' specified)\n"
"      dosagesum: Sum of named allele dosages.\n"
"      scoreavgs: Score averages.\n"
"      scoresums: Score sums.\n"
"    The default is maybesid,phenos,nmissallele,dosagesum,scoreavgs.\n\n"
	       );
    help_print("genotyping-rate", &help_ctrl, 1,
"  --genotyping-rate <dosage>\n"
"    Report genotyping rate in log (this was automatic in PLINK 1.x).\n\n"
	       );
    help_print("validate", &help_ctrl, 1,
"  --validate\n"
"    Validates all variant records in a .pgen file.\n\n"
	       );
    help_print("zst-decompress", &help_ctrl, 1,
"  --zst-decompress [.zst file] {output filename}\n"
"    Decompress a Zstd-compressed file.  If no output filename is specified, the\n"
"    file is decompressed to standard output.\n"
"    This cannot be used with any other flags, and does not cause a log file to\n"
"    be generated.\n\n"
	       );
    if (!param_ct) {
      fputs(
"The following other flags are supported.\n"
// tbd: document order of operations
, stdout);
    }
    help_print("script\trerun", &help_ctrl, 0,
"  --script [fname]   : Include command-line options from file.\n"
"  --rerun {log}      : Rerun commands in log (default '" PROG_NAME_STR ".log').\n"
	       );
    help_print("version", &help_ctrl, 0,
"  --version          : Display only version number before exiting.\n"
	       );
    help_print("silent", &help_ctrl, 0,
"  --silent           : Suppress output to console.\n"
	       );
    help_print("input-missing-genotype\tmissing-genotype", &help_ctrl, 0,
"  --input-missing-genotype [c] : '.' is always interpreted as a missing\n"
"                                 genotype code in input files.  By default, '0'\n"
"                                 also is; you can change this second missing\n"
"                                 code with --input-missing-genotype.\n"
	       );
    help_print("vcf\tbcf\tbgen\tdouble-id\tconst-fid\tid-delim", &help_ctrl, 0,
"  --double-id        : Set both FIDs and IIDs to the VCF/.bgen sample ID.\n"
"  --const-fid {ID}   : Set all FIDs to the given constant (default '0').\n"
"  --id-delim {d}     : Parse sample IDs as [FID][d][IID] (or\n"
"                       [FID][d][IID][d][SID] when delimiter appears twice).\n"
"                       Default delimiter is '_'.\n"
	       );
    help_print("idspace-to\tvcf\tbcf\tbgen\tid-delim\tvcf-idspace-to", &help_ctrl, 0,
"  --idspace-to [c]   : Convert spaces in VCF/.bgen sample IDs to the given\n"
"                       character.\n"
	       );
    help_print("vcf\tbcf\tvcf-half-call\tvcf-min-gq\tvcf-min-dp\tvcf-require-gt", &help_ctrl, 0,
"  --vcf-require-gt   : Skip variants with no GT field.\n"
"  --vcf-min-gq [val] : No-call genotypes when GQ is present and below the\n"
"                       threshold.\n"
"  --vcf-min-dp [val] : No-call genotypes when DP is present and below the\n"
"                       threshold.\n"
"  --vcf-half-call [] : Specify how '0/.' and similar VCF GT values should be\n"
"                       handled.  The following four modes are supported:\n"
"                       * 'error'/'e' (default) errors out and reports line #.\n"
"                       * 'haploid'/'h' treats them as haploid calls.\n"
"                       * 'missing'/'m' treats them as missing.\n"
"                       * 'reference'/'r' treats the missing value as 0.\n"
	       );
    help_print("oxford-single-chr\tdata\tgen", &help_ctrl, 0,
"  --oxford-single-chr [chr name]  : Specify single-chromosome .gen file with\n"
"                                    ignorable first column.\n"
	       );
    // any need to keep --hard-call-threshold random?  postpone it for now...
    help_print("hard-call-threshold\tgen\tbgen\tdata\timport-dosage", &help_ctrl, 0,
"  --hard-call-threshold [val]     : When importing dosage data, a hardcall is\n"
"                                    normally saved when the distance from the\n"
"                                    nearest hardcall, defined as\n"
"                                      0.5 * sum_i |x_i - round(x_i)|\n"
"                                    (where the x_i's are 0..2 allele dosages),\n"
"                                    is not greater than 0.1.  You can adjust\n"
"                                    this threshold by providing a numeric\n"
"                                    parameter to --hard-call-threshold.\n"
"                                    You can also use this with --make-{b}pgen\n"
"                                    to alter the saved hardcalls while leaving\n"
"                                    the dosages untouched.\n"
	       );
    help_print("dosage-erase-threshold\timport-dosage-certainty\tgen\tbgen\tdata\tvcf\tbcf\timport-dosage", &help_ctrl, 0,
"  --dosage-erase-threshold [val]  : --hard-call-threshold normally preserves\n"
"                                    the original dosages, and several PLINK 2.x\n"
"                                    commands use them when they're available.\n"
"                                    Use --dosage-erase-threshold to make PLINK\n"
"                                    erase dosages and keep only hardcalls when\n"
"                                    distance-from-hardcall <= the given level.\n"
"  --import-dosage-certainty [val] : The PLINK 2.0 file format currently\n"
"                                    supports a single dosage for each allele.\n"
"                                    Some other dosage file formats include a\n"
"                                    separate probability for every possible\n"
"                                    genotype, e.g. {P(0/0)=0.2, P(0/1)=0.52,\n"
"                                    P(1/1)=0.28}, a highly uncertain call that\n"
"                                    is nevertheless treated as a hardcall under\n"
"                                    '--hard-call-threshold 0.1'.  To make PLINK\n"
"                                    treat a dosage as missing whenever the\n"
"                                    largest probability is less than a\n"
"                                    threshold, use --import-dosage-certainty.\n"
	       );
    help_print("missing-code\tmissing_code\tdata\tsample", &help_ctrl, 0,
"  --missing-code {string list}    : Comma-delimited list of missing phenotype\n"
"    (alias: --missing_code)         values for Oxford-format import (default\n"
"                                    'NA').\n"
	       );
    help_print("allow-extra-chr\taec", &help_ctrl, 0,
"  --allow-extra-chr  : Permit unrecognized chromosome codes (alias --aec).\n"
	       );
    // possible todo: nonhuman PARs?
    help_print("chr-set\tcow\tdog\thorse\thound\tmouse\trice\tsheep\tautosome-num\thuman\tchr-override", &help_ctrl, 0,
"  --chr-set [autosome ct] <no-x> <no-y> <no-xy> <no-mt> :\n"
"    Specify a nonhuman chromosome set.  The first parameter sets the number of\n"
"    diploid autosome pairs if positive, or haploid chromosomes if negative.\n"
"    Given diploid autosomes, the remaining modifiers indicate the absence of\n"
"    the named non-autosomal chromosomes.\n"
"  --cow/--dog/--horse/--mouse/--rice/--sheep : Shortcuts for those species.\n"
"  --autosome-num [val]  : Alias for '--chr-set [value] no-y no-xy no-mt'.\n"
"  --human               : Explicitly specify human chromosome set, and make\n"
"                          output .pvar/VCF files include a ##chrSet header\n"
"                          line.  (.pvar/VCF output files automatically include\n"
"                          ##chrSet when a nonhuman set is specified.)\n"
"  --chr-override <file> : By default, if --chr-set/--autosome-num/--human/etc.\n"
"                          conflict with an input file ##chrSet header line,\n"
"                          PLINK will error out.  --chr-override with no\n"
"                          parameter causes the command line to take precedence;\n"
"                          '--chr-override file' defers to the file.\n"
	       );
    help_print("biallelic-only\tvar-min-qual\tvar-filter\tvcf-min-qual\tvcf-filter", &help_ctrl, 0,
"  --biallelic-only <strict> <list> : Skip variants with 2+ alt. alleles.\n"
"  --var-min-qual [val]             : Skip variants with low/missing QUAL.\n"
"  --var-filter {exception(s)...}   : Skip variants which have FILTER failures.\n"
	       );
    /*
    help_print("allow-no-samples\tallow-no-vars", &help_ctrl, 0,
"  --allow-no-samples : Allow the input fileset to contain no samples.\n"
"  --allow-no-vars    : Allow the input fileset to contain no variants.\n"
	       );
    */
    help_print("pheno\tpheno-name", &help_ctrl, 0,
"  --pheno [filename] : Specify additional phenotype/covariate file.\n"
"  --pheno-name [...] : Only load the designated phenotype(s) from the\n"
"                       --pheno (if one was specified) or .psam (if no --pheno)\n"
"                       file.  Separate multiple names with spaces or commas,\n"
"                       and use dashes to designate ranges.\n"
	       );
    help_print("input-missing-phenotype\t1\tmissing-catname\tmissing-phenotype", &help_ctrl, 0,
"  --input-missing-phenotype [v] : Set number to treat as a missing phenotype in\n"
"                                  input files (default -9).\n"
"  --1                           : Expect case/control phenotypes in input files\n"
"                                  to be coded as 0 = control, 1 = case, instead\n"
"                                  of the usual 0 = missing, 1 = ctrl, 2 = case.\n"
"  --missing-catname [str]       : Set missing-categorical-phenotype string\n"
"                                  (case-sensitive, default 'NONE').\n"
	       );
    help_print("covar\tcovar-name", &help_ctrl, 0,
"  --covar [filename] : Specify additional covariate file.\n"
"  --covar-name [...] : Only load the designated covariate(s) from the\n"
"                       --covar (if one was specified), --pheno (if no --covar),\n"
"                       or .psam (if no --covar or --pheno) file.\n"
	       );
    help_print("within\tmwithin\tfamily\tfamily-missing-catname", &help_ctrl, 0,
"  --within [f] {new pheno name} : Import a PLINK 1.x categorical phenotype.\n"
"                                  (Phenotype name defaults to 'CATPHENO'.)\n"
"                                  * If any numeric values are present, ALL\n"
"                                    values must be numeric.  In that case, 'C'\n"
"                                    is added in front of all category names.\n"
"                                  * 'NA' is treated as a missing value.\n"
"  --mwithin [n]                 : Load --within categories from column n+2.\n"
"  --family {new pheno name}     : Create a categorical phenotype from FID.\n"
"                                  Restrictions on and handling of numeric\n"
"                                  values are the same as for --within.\n"
"  --family-missing-catname [nm] : Make --family treat the specified FID as\n"
"                                  missing.\n"
	       );
    help_print("keep\tremove\tkeep-fam\tremove-fam", &help_ctrl, 0,
"  --keep <sid> [fn...]  : Exclude all samples not named in a file.\n"
"  --remove <sid> [f...] : Exclude all samples named in a file.\n"
"  --keep-fam [fname...] : Exclude all families not named in a file.\n"
"  --remove-fam [fn...]  : Exclude all families named in a file.\n"
	       );
    help_print("extract\texclude\trange", &help_ctrl, 0,
"  --extract <range> [f...] : Exclude all variants not named in a file.\n"
"  --exclude <range> [f...] : Exclude all variants named in a file.\n"
	       );
    help_print("keep-cats\tkeep-cat-names\tkeep-cat-pheno\tremove-cats\tremove-cat-names\tremove-cat-pheno\tkeep-clusters\tkeep-cluster-names\tremove-clusters\tremove-cluster-names", &help_ctrl, 0,
"  --keep-cats [filename]   : These can be used individually or in combination\n"
"  --keep-cat-names [nm...]   to define a list of categories to keep; all\n"
"                             samples not in one of the named categories are\n"
"                             excluded.  Use spaces to separate category names\n"
"                             for --keep-cat-names.  Use the --missing-catname\n"
"                             value (default 'NONE') to refer to the group of\n"
"                             uncategorized samples.\n"
"  --keep-cat-pheno [pheno] : If more than one categorical phenotype is loaded,\n"
"                             or you wish to filter on a categorical covariate,\n"
"                             --keep-cat-pheno must be used to specify which\n"
"                             phenotype/covariate --keep-cats and\n"
"                             --keep-cat-names apply to.\n"
"  --remove-cats [filename] : Exclude all categories named in the file.\n"
"  --remove-cat-names [...] : Exclude named categories.\n"
"  --remove-cat-pheno [phe] : Specify pheno for --remove-cats/remove-cat-names.\n"
	       );
    help_print("split-cat-pheno\tdummy-coding\tloop-assoc", &help_ctrl, 0,
"  --split-cat-pheno <omit-last> <covar-01> {cat. pheno/covar name(s)...} :\n"
"    Split n-category phenotype(s) into n (or n-1, with 'omit-last') binary\n"
"    phenotypes, with names of the form [orig. pheno name]=[category name].  (As\n"
"    a consequence, affected phenotypes and categories are not permitted to\n"
"    contain the '=' character.)\n"
"    * This happens after all sample filters.\n"
"    * If no phenotype or covariate names are provided, all categorical\n"
"      phenotypes (but not covariates) are processed.\n"
"    * By default, generated covariates are coded as 1=false, 2=true.  To code\n"
"      them as 0=false, 1=true instead, add the 'covar-01' modifier.\n"
	       );
    help_print("variance-standardize\tcovar-variance-standardize\tquantile-normalize\tpheno-quantile-normalize\tcovar-quantile-normalize\tstandard-beta\tglm\tlinear\tlogistic", &help_ctrl, 0,
"  --variance-standardize {pheno/covar name(s)...}\n"
"  --covar-variance-standardize {covar name(s)...} :\n"
"    Linearly transform named covariates (and quantitative phenotypes, if\n"
"    --variance-standardize) to mean-zero, variance 1.  If no parameters are\n"
"    provided, all possible phenotypes/covariates are affected.\n"
"  --quantile-normalize {...}       : Force named covariates and quantitative\n"
"  --pheno-quantile-normalize {...}   phenotypes to a N(0,1) distribution,\n"
"  --covar-quantile-normalize {...}   preserving only the original rank orders.\n"
	       );
    help_print("chr\tnot-chr", &help_ctrl, 0,
"  --chr [chr(s)...]  : Exclude all variants not on the given chromosome(s).\n"
"                       Valid choices for humans are 0 (unplaced), 1-22, X, Y,\n"
"                       XY, MT, PAR1, and PAR2.  Separate multiple chromosomes\n"
"                       with spaces and/or commas, and use a dash (no adjacent\n"
"                       spaces permitted) to denote a range, e.g.\n"
"                       '--chr 1-4, 22, par1, x, par2'.\n"
"  --not-chr [...]    : Reverse of --chr (exclude variants on listed\n"
"                       chromosomes).\n"
	       );
    help_print("autosome\tautosome-par\tautosome-xy\tchr\tnot-chr", &help_ctrl, 0,
"  --autosome         : Exclude all non-autosomal variants.\n"
"  --autosome-par     : Exclude all non-autosomal variants, except those in a\n"
"                       pseudo-autosomal region.\n"
	       );
    help_print("snps-only", &help_ctrl, 0,
"  --snps-only <just-acgt> : Exclude non-SNP variants.  By default, SNP = all\n"
"                            allele codes are single-character; 'just-acgt'\n"
"                            restricts SNP codes to {A,C,G,T,a,c,g,t,[missing]}.\n"
	       );
    // best to only support --chr with --from-bp/--to-bp/etc., now that
    // finalize_chrset() is deferred
    help_print("from\tto\tsnp\twindow\tfrom-bp\tto-bp\tfrom-kb\tto-kb\tfrom-mb\tto-mb\texclude-snp\textract-snp", &help_ctrl, 0,
"  --from [var ID]    : Use ID(s) to specify a variant range to load.  When used\n"
"  --to   [var ID]      together, both variants must be on the same chromosome.\n"
"                       (--snps can be used to specify intervals which cross\n"
"                       chromosome boundaries.)\n"
"  --snp  [var ID]    : Specify a single variant to load.\n"
"  --exclude-snp [ID] : Specify a single variant to exclude.\n"
"  --window  [kbs]    : With --snp/--exclude-snp, loads/excludes all variants\n"
"                       within half the specified kb distance of the named one.\n"
"  --from-bp [pos]    : Use base-pair coordinates to define a variant range to\n"
"  --to-bp   [pos]      load.\n"
"  --from-kb [pos]      * You must use these with --chr, specifying a single\n"
"  --to-kb   [pos]        chromosome.\n"
"  --from-mb [pos]      * Decimals and negative numbers are permitted.\n"
"  --to-mb   [pos]      * The --to-bp(/-kb/-mb) position is no longer permitted\n"
"                         to be smaller than the --from-bp position.\n"
	       );
    help_print("snps\texclude-snps", &help_ctrl, 0,
"  --snps [var IDs...]  : Use IDs to specify variant range(s) to load or\n"
"  --exclude-snps [...]   exclude.  E.g. '--snps rs1111-rs2222, rs3333, rs4444'.\n"
	       );
    help_print("force-intersect\textract\tfrom\tto\tsnp\tfrom-bp\tto-bp\tfrom-kb\tto-kb\tfrom-mb\tto-mb\textract-snp\tsnps", &help_ctrl, 0,
"  --force-intersect    : PLINK 2 normally errors out when multiple variant\n"
"                         inclusion filters (--extract, --from/--to,\n"
"                         --from-bp/--to-bp, --snp, --snps) are specified.\n"
"                         --force-intersect allows the run to proceed; the set\n"
"                         intersection will be taken.\n"
	       );
    help_print("geno\tmind\toblig-clusters\toblig-missing", &help_ctrl, 0,
"  --geno {val} <dosage | hh-missing>\n"
"  --mind {val} <dosage | hh-missing> : \n"
"    Exclude variants (--geno) and/or samples (--mind) with missing call\n"
"    frequencies greater than a threshold (default 0.1).  (Note that the default\n"
"    threshold is only applied if --geno/--mind is invoked without a parameter;\n"
"    when --geno/--mind is not invoked, no missing call frequency ceiling is\n""    enforced at all.  Other inclusion/exclusion default thresholds work the\n"
"    same way.)\n"
"    By default, when a dosage is present but a hardcall is not, the genotype is\n"
"    treated as missing; add the 'dosage' modifier to treat this case as\n"
"    nonmissing.  Alternatively, you can use 'hh-missing' to also treat\n"
"    heterozygous haploid calls as missing.\n"
	       );
    /*
    help_print("oblig-clusters\toblig-missing", &help_ctrl, 0,
"  --oblig-missing [f1] [f2] : Specify blocks of missing genotype calls for\n"
"                              --geno/--mind to ignore.  The first file should\n"
"                              have variant IDs in the first column and block\n"
"                              IDs in the second, while the second file should\n"
"                              have FIDs in the first column, IIDs in the\n"
"                              second, and block IDs in the third.\n"
	       );
    */
    help_print("require-pheno\trequire-covar\tprune", &help_ctrl, 0,
"  --require-pheno {name(s)...} : Remove samples missing any of the named\n"
"  --require-covar {name(s)...}   phenotype(s)/covariate(s).  If no parameters\n"
"                                 are provided, all phenotype(s)/covariate(s)\n"
"                                 must be present.\n"
	       );
    help_print("maf\tmax-maf\tmac\tmin-ac\tmax-mac\tmax-ac", &help_ctrl, 0,
"  --maf {freq}       : Exclude variants with nonmajor allele frequency lower\n"
"                       than a threshold (default 0.01).\n"
"  --max-maf [freq]   : Exclude variants with MAF greater than the threshold.\n"
"  --mac [ct]         : Exclude variants with nonmajor allele dosage lower than\n"
"                       the given threshold.\n"
"  --max-mac [ct]     : Exclude variants with nonmajor allele dosage greater than\n"
"                       the given threshold.\n"
	       );
    help_print("maf-succ", &help_ctrl, 0,
"  --maf-succ         : Rule of succession allele frequency estimation (used in\n"
"                       EIGENSOFT).  Given a j observations of one allele and k\n"
"                       observations of the other for a biallelic variant, infer\n"
"                       allele frequencies of (j+1) / (j+k+2) and\n"
"                       (k+1) / (j+k+2), rather than the default j / (j+k) and\n"
"                       k / (j+k).\n"
"                       Note that this does not affect --freq's output.\n"
	       );
    help_print("read-freq", &help_ctrl, 0,
"  --read-freq [file] : Load allele frequency estimates from the given --freq or\n"
"                       --geno-counts (or PLINK 1.9 --freqx) report, instead of\n"
"                       imputing them from the immediate dataset.\n"
	       );
// todo: something like <check-ctrls>/<check-ctrl=[case/ctrl phenotype name]>
// and maybe <ctrls-only>/<ctrl-only=[case/ctrl phenotype name]>
    help_print("hwe\tmach-r2-filter", &help_ctrl, 0,
"  --hwe [p] <midp> <keep-fewhet> : Exclude variants with Hardy-Weinberg\n"
"                                   equilibrium exact test p-values below a\n"
"                                   threshold.\n"
"                                   * By default, only founders are considered.\n"
"                                   * chrX p-values are now computed using\n"
"                                     Graffelman and Weir's method.\n"
"                                   * With 'keep-fewhet', variants which fail\n"
"                                     the test in the too-few-hets direction are\n"
"                                     not excluded.  (On chrX, this uses the\n"
"                                     ratio between the Graffelman/Weir p-value\n"
"                                     and the female-only p-value.)\n"
"                                   * There is currently no special handling of\n"
"                                     case/control phenotypes.\n"
"  --mach-r2-filter {min} {max}   : Exclude variants with MaCH\n"
"                                   empirical-theoretical variance ratio outside\n"
"                                   of [min, max] (defaults 0.1 and 2.0).\n"
"                                   * For multiallelic variants, only the\n"
"                                     ref-nonref dimension is considered.\n"
"                                   * If a single parameter is provided, it is\t"
"                                     treated as the minimum.\n"
	       );
    help_print("keep-females\tkeep-males\tkeep-nosex\tremove-females\tremove-males\tremove-nosex\tfilter-males\tfilter-females", &help_ctrl, 0,
"  --keep-females     : Exclude male and unknown-sex samples.\n"
"  --keep-males       : Exclude female and unknown-sex samples.\n"
"  --keep-nosex       : Exclude all known-sex samples.\n"
"  --remove-females   : Exclude female samples.\n"
"  --remove-males     : Exclude male samples.\n"
"  --remove-nosex     : Exclude unknown-sex samples.\n"
	       );
    help_print("keep-founders\tkeep-nonfounders\tfilter-founders\tfilter-nonfounders\tgeno-counts", &help_ctrl, 0,
"  --keep-founders    : Exclude nonfounder samples.\n"
"  --keep-nonfounders : Exclude founder samples.\n"
	       );
    // possible todo: allow or/and of multiple predicates
    // best if syntax allows for '=' character inside phenotype/covariate
    // names, though...
    help_print("keep-if\tremove-if\tfilter-cases\tfilter-controls", &help_ctrl, 0,
"  --keep-if [pheno/covar] [op] [val] : Exclude samples which don't/do satisfy a\n"
"  --remove-if [pheno/covar] [op] [v]   comparison predicate, e.g.\n"
"                                         --keep-if PHENO1 == case\n"
"                                       Unless the operator is !=, the predicate\n"
"                                       always evaluates to false when the\n"
"                                       phenotype/covariate is missing.\n"
	       );
    help_print("nonfounders\tfreq\thardy\thwe", &help_ctrl, 0,
"  --nonfounders      : Include nonfounders in allele freq/HWE calculations.\n"
	       );
    help_print("output-chr", &help_ctrl, 0,
"  --output-chr [MT code] : Set chromosome coding scheme in output files by\n"
"                           providing the desired human mitochondrial code.\n"
"                           Options are '26', 'M', 'MT', '0M', 'chr26', 'chrM',\n"
"                           and 'chrMT'; default is now 'MT' (note that this is\n"
"                           a change from PLINK 1.x, which defaulted to '26').\n"
	       );
    help_print("output-missing-genotype\toutput-missing-phenotype\tmissing-genotype\tmissing-phenotype", &help_ctrl, 0,
"  --output-missing-genotype [ch] : Set the code used to represent missing\n"
"                                   genotypes in output files (default '.').\n"
"  --output-missing-phenotype [s] : Set the string used to represent missing\n"
"                                   phenotypes in output files (default 'NA').\n"
	       );
    /*
    help_print("sort-vars", &help_ctrl, 0,
"  --sort-vars {mode}      : Sort variants by chromosome, then position, then\n"
"                            ID.  The following string orders are supported:\n"
"                            * 'natural'/'n': Natural sort (default).\n"
"                            * 'ascii'/'a': ASCII.\n"
"                            This must be used with --make-{b}pgen/--make-bed.\n"
	       );
    */
    help_print("set-hh-missing\tset-mixed-mt-missing", &help_ctrl, 0,
"  --set-hh-missing        : Make --make-{b}pgen/--make-bed set heterozygous\n"
"                            haploid and female chrY genotypes to missing.\n"
"                            (Unlike PLINK 1.x, this does not change unknown-sex\n"
"                            chrY genotypes.)\n"
"  --set-mixed-mt-missing  : Make --make-{b}pgen/--make-bed set mixed MT\n"
"                            genotypes to missing.\n"
	       );
    help_print("split-par\tmerge-par\tsplit-x\tmerge-x", &help_ctrl, 0,
"  --split-par [bp1] [bp2] : Changes chromosome code of all X chromosome\n"
"  --split-par [build]       variants with bp position <= bp1 to PAR1, and those\n"
"                            with position >= bp2 to PAR2.  The following build\n"
"                            codes are supported as shorthand:\n"
"                            * 'b36'/'hg18' = NCBI 36, 2709521/154584237\n"
"                            * 'b37'/'hg19' = GRCh37, 2699520/154931044\n"
"                            * 'b38'/'hg38' = GRCh38, 2781479/155701383\n"
"  --merge-par             : Merge PAR1/PAR2 back with X.  Requires PAR1 to be\n"
"                            positioned immediately before X, and PAR2 to be\n"
"                            immediately after X.  (Should *not* be used with\n"
"                            \"--export vcf\", since it causes male\n"
"                            homozygous/missing calls in PAR1/PAR2 to be\n"
"                            reported as haploid.)\n"
	       );
    help_print("set-all-var-ids\tset-missing-var-ids\tnew-id-max-allele-len\tmissing-var-code", &help_ctrl, 0,
"  --set-missing-var-ids [t]  : Given a template string with a '@' where the\n"
"  --set-all-var-ids [t]        chromosome code should go and '#' where the bp\n"
"                               coordinate belongs, --set-missing-var-ids\n"
"                               assigns chromosome-and-bp-based IDs to unnamed\n"
"                               variants, while --set-all-var-ids resets all\n"
"                               IDs.\n"
"                               You may also use '$r'/'$a' to refer to the\n"
"                               ref and alt1 alleles, or '$1'/'$2' to refer to\n"
"                               them in alphabetical order.\n"
"  --new-id-max-allele-len [len] <error | missing | truncate> :\n"
"    Specify maximum number of leading characters from allele codes to include\n"
"    in new variant IDs, and behavior on longer codes (defaults 23, error).\n"
"  --missing-var-code [str]   : Change unnamed variant code for\n"
"                               --set-[missing/all]-var-ids (default '.').\n"
	       );
    help_print("update-sex", &help_ctrl, 0,
"  --update-sex [f] {n} : Update sexes.  Sex (1/M/m = male, 2/F/f = female, 0 =\n"
"                         missing) is loaded from column n+2 (default n is 1).\n"
	       );
    // don't make --real-ref-alleles apply to e.g. Oxford import, since
    // explicit 'ref-first'/'ref-second' modifiers are clearer
    help_print("real-ref-alleles\tmaj-ref\tkeep-allele-order", &help_ctrl, 0,
"  --real-ref-alleles : Treat A2 alleles in a PLINK 1.x fileset as actual ref\n"
"                       alleles; otherwise they're flagged as provisional.\n"
"  --maj-ref <force>  : Set major alleles to reference, like PLINK 1.x\n"
"                       automatically did.  (Note that this is now opt-in rather\n"
"                       than opt-out; --keep-allele-order is no longer necessary\n"
"                       to prevent allele-swapping.)  By default, this only\n"
"                       affects variants with \"provisional reference\" flags;\n"
"                       add 'force' to override validated reference alleles as\n"
"                       well.\n"
"                       All new reference alleles are marked as provisional.\n"
	       );
    help_print("indiv-sort", &help_ctrl, 0,
"  --indiv-sort [m] <sid> {f} : Specify FID/IID(/SID) sort order for merge and\n"
"                               --make-{b}pgen/--make-bed.  The following four\n"
"                               modes are supported:\n"
"                               * 'none'/'0' keeps samples in the order they\n"
"                                 were loaded.  Default for non-merge.\n"
"                               * 'natural'/'n' invokes \"natural sort\", e.g.\n"
"                                 'id2' < 'ID3' < 'id10'.  Default when merging.\n"
"                               * 'ascii'/'a' sorts in ASCII order, e.g.\n"
"                                 'ID3' < 'id10' < 'id2'.\n"
"                               * 'file'/'f' uses the order in the given file\n"
"                                 (named in the last parameter).  The 'sid'\n"
"                                 modifier has the usual effect when this mode\n"
"                                 is requested.\n"
	       );
    help_print("make-king\tmake-king-table\tking-table-filter", &help_ctrl, 0,
"  --king-table-filter [min]  : Specify minimum kinship coefficient for\n"
"                               inclusion in --make-king-table report.\n"
	       );
    help_print("glm\tlinear\tlogistic\tcondition\tcondition-list\tparameters\ttests", &help_ctrl, 0,
"  --condition [var ID] <dominant | recessive> : Add one variant's alt1 dosages\n"
"                                                as a --glm covariate.\n"
"  --condition-list [f] <dominant | recessive> : Add all variants in the file as\n"
"                                                --glm covariates.\n"
"  --parameters [...] : Include only the given covariates/interactions in the\n"
"                       --glm model, identified by a list of 1-based indices\n"
"                       and/or ranges of them.\n"
	       /*
"  --tests [...]      : Perform a (joint) test on the specified term(s) in the\n"
"  --tests all          --glm model, identified by 1-based indices and/or ranges\n"
"                       of them.  If permutation was requested, it is based on\n"
"                       this test.\n"
"                       * Note that, when --parameters is also present, the\n"
"                         indices refer to the terms remaining AFTER pruning by\n"
"                         --parameters.\n"
"                       * You can use '--tests all' to include all terms.\n"
	       */
	       );
    help_print("glm\tlinear\tlogistic\tvif\tmax-corr", &help_ctrl, 0,
"  --vif [max VIF]    : Set VIF threshold for --glm multicollinearity check\n"
"                       (default 50).  For case/control phenotypes, only the\n"
"                       covariates are checked, and the entire phenotype is\n"
"                       skipped if a VIF is too high.\n"
"  --max-corr [val]   : Skip --glm regression when the absolute value of the\n"
"                       correlation between two predictors exceeds this value\n"
"                       (default 0.999).  For case/control phenotypes, only\n"
"                       covariates are checked.\n"
	       );
    help_print("glm\tlinear\tlogistic\tscore\txchr-model", &help_ctrl, 0,
"  --xchr-model [m]   : Set the chrX --glm/--score model.\n"
"                       * '0' = skip chrX.\n"
"                       * '1' = add sex as a covar on chrX, code males 0..1.\n"
"                       * '2' (default) = chrX sex covar, code males 0..2.\n"
"                       (Use the --glm 'interaction' modifier to test for\n"
"                       interaction between genotype and sex.)\n"
	       );
    /*
    help_print("adjust", &help_ctrl, 0,
"  --adjust <gc> <log10> <cols=[column set descriptor]> :\n"
"    For each association test, report some multiple-testing corrections, sorted\n"
"    in increasing-p-value order.\n"
"    * 'gc' causes genomic-controlled p-values to be used in the formulas.\n"
"    * 'log10' causes negative base 10 logs of p-values to be reported, instead\n"
"      of raw p-values.\n"
"    The following column sets are supported:\n"
"      chrom: Chromosome ID.\n"
"      pos: Base-pair coordinate.\n"
"      (ID is always present, and positioned here.)\n"
"      ref: Reference allele.\n"
"      alt1: Alternate allele 1.\n"
"      alt: All alternate alleles, comma-separated.\n"
"      unadj: Unadjusted p-value.\n"
"      gc: Devlin & Roeder (1999) genomic control corrected p-value (additive\n"
"          models only).\n"
"      qq: P-value quantile.\n"
"      bonf: Bonferroni correction.\n"
"      holm: Holm-Bonferroni (1979) adjusted p-value.\n"
"      sidakss: Sidak single-step adjusted p-value.\n"
"      sidaksd: Sidak step-down adjusted p-value.\n"
"      fdrbh: Benjamini & Hochberg (1995) step-up false discovery control.\n"
"      fdrby: Benjamini & Yekutieli (2001) step-up false discovery control.\n"
"    Default set is chrom,unadj,gc,bonf,holm,sidakss,sidaksd,fdrbh,fdrby.\n"
	       );
    help_print("adjust\tlambda", &help_ctrl, 0,
"  --lambda           : Set genomic control lambda for --adjust.\n"
	       );
    */
    help_print("ci\tlinear\tlogistic", &help_ctrl, 0,
"  --ci [size]        : Report confidence ratios for odds ratios/betas.\n"
	       );
    help_print("pfilter", &help_ctrl, 0,
"  --pfilter [val]    : Filter out assoc. test results with higher p-values.\n"
	       );
    /*
    help_print("aperm", &help_ctrl, 0,
"  --aperm [min perms - 1] {max perms} {alpha} {beta} {init interval} {slope} :\n"
"    Set up to six parameters controlling adaptive permutation tests.\n"
"    * The first two control the minimum and maximum number of permutations that\n"
"      may be run for each variant; default values are 5 and 1000000.\n"
"    * The next two control the early termination condition.  A\n"
"      100% * (1 - beta/2T) confidence interval is calculated for each empirical\n"
"      p-value, where T is the total number of variants; whenever this\n"
"      confidence interval doesn't contain alpha, the variant is exempted from\n"
"      further permutation testing.  Default values are 0 and 1e-4.\n"
"    * The last two control when the early termination condition is checked.  If\n"
"      a check occurs at permutation #p, the next check occurs after\n"
"      [slope]p + [init interval] more permutations (rounded down).  Default\n"
"      initial interval is 1, and default slope is 0.001.\n"
	       );
    help_print("mperm-save\tmperm-save-all", &help_ctrl, 0,
"  --mperm-save       : Save best max(T) permutation test statistics.\n"
"  --mperm-save-all   : Save all max(T) permutation test statistics.\n"
	       );
    */
    help_print("score-col-nums\tscore", &help_ctrl, 0,
"  --score-col-nums [...] : Process all the specified coefficient columns in the\n"
"                           --score file, identified by 1-based indexes and/or\n"
"                           ranges of them.\n"
	       );
    help_print("parallel", &help_ctrl, 0,
"  --parallel [k] [n] : Divide the output matrix into n pieces, and only compute\n"
"                       the kth piece.  The primary output file will have the\n"
"                       piece number included in its name, e.g. plink2.king.13\n"
"                       or plink2.king.13.zst if k is 13.  Concatenating these\n"
"                       files in order will yield the full matrix of interest.\n"
"                       (Yes, this can be done before decompression.)\n"
"                       N.B. This generally cannot be used to directly write a\n"
"                       symmetric square matrix.  Choose square0 or triangle\n"
"                       shape instead, and postprocess as necessary.\n"
	       );
    help_print("memory\tseed", &help_ctrl, 0,
"  --memory [val] <require> : Set size, in MB, of initial workspace malloc\n"
"                             attempt.  To error out instead of reducing the\n"
"                             request size when the initial attempt fails, add\n"
"                             the 'require' modifier.\n"
	       );
    help_print("threads\tnum_threads\tthread-num\tseed", &help_ctrl, 0,
"  --threads [val]    : Set maximum number of compute threads.\n"
	       );
    help_print("seed", &help_ctrl, 0,
"  --seed [val...]    : Set random number seed(s).  Each value must be an\n"
"                       integer between 0 and 4294967295 inclusive.\n"
"                       Note that --threads and \"--memory require\" may also be\n"
"                       needed to reproduce some randomized runs.\n"
	       );
    help_print("output-min-p", &help_ctrl, 0,
"  --output-min-p [p] : Specify minimum p-value to write to reports.\n"
	       );
    help_print("debug", &help_ctrl, 0,
"  --debug            : Use slower, more crash-resistant logging method.\n"
	       );
    help_print("warning-errcode", &help_ctrl, 0,
"  --warning-errcode  : Return a nonzero error code to the OS when a run\n"
"                       completes with warning(s).\n"
	       );
    if (!param_ct) {
      fputs(
"\nPrimary methods paper:\n"
"Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ (2015)\n"
"Second-generation PLINK: rising to the challenge of larger and richer datasets.\n"
"GigaScience, 4.\n"
, stdout);
    }
  } while (help_ctrl.iters_left--);
  if (help_ctrl.unmatched_ct) {
    net_unmatched_ct = help_ctrl.unmatched_ct;
    printf("\nNo help entr%s for", (help_ctrl.unmatched_ct == 1)? "y" : "ies");
    col_num = (help_ctrl.unmatched_ct == 1)? 17 : 19;
    arg_uidx = 0;
    // er, should replace the \n logic with a wordwrap() call
    while (help_ctrl.unmatched_ct) {
      arg_uidx = next_unset_unsafe(help_ctrl.all_match_arr, arg_uidx);
      help_ctrl.unmatched_ct--;
      if (help_ctrl.unmatched_ct) {
	if (net_unmatched_ct == 2) {
	  if (help_ctrl.param_slens[arg_uidx] + col_num > 76) {
	    putc_unlocked('\n', stdout);
	    col_num = 2 + help_ctrl.param_slens[arg_uidx];
	  } else {
	    putc_unlocked(' ', stdout);
	    col_num += 3 + help_ctrl.param_slens[arg_uidx];
	  }
	  putc_unlocked('\'', stdout);
	  fputs(argv[arg_uidx], stdout);
	  putc_unlocked('\'', stdout);
	} else {
	  if (help_ctrl.param_slens[arg_uidx] + col_num > 75) {
	    putc_unlocked('\n', stdout);
	    col_num = 3 + help_ctrl.param_slens[arg_uidx];
	  } else {
	    putc_unlocked(' ', stdout);
	    col_num += 4 + help_ctrl.param_slens[arg_uidx];
	  }
	  putc_unlocked('\'', stdout);
	  fputs(argv[arg_uidx], stdout);
          fputs("',", stdout);
	}
	if (help_ctrl.unmatched_ct == 1) {
	  if (col_num > 76) {
	    fputs("\nor", stdout);
	    col_num = 2;
	  } else {
	    fputs(" or", stdout);
	    col_num += 3;
	  }
	}
      } else {
	putc_unlocked((help_ctrl.param_slens[arg_uidx] + col_num > 75)? '\n' : ' ', stdout);
	putc_unlocked('\'', stdout);
        fputs(argv[arg_uidx], stdout);
        fputs("\'.\n", stdout);
      }
      arg_uidx++;
    }
  }
  if (param_ct) {
    while (0) {
    disp_help_ret_NOMEM:
      reterr = kPglRetNomem;
    }
    free_cond(help_ctrl.param_slens);
    free_cond(help_ctrl.all_match_arr);
    if (help_ctrl.argv && (help_ctrl.argv != argv)) {
      free(help_ctrl.argv);
    }
  }
  return reterr;
}

#ifdef __cplusplus
}
#endif

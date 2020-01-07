// This file is part of PLINK 1.90, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "plink_common.h"

uint32_t edit1_match(uint32_t len1, char* s1, uint32_t len2, char* s2) {
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
	diff_found++;
      }
      pos++;
    }
  } else if (len1 == len2 - 1) {
    do {
      if (s1[pos - diff_found] != s2[pos]) {
	if (diff_found) {
	  return 0;
	}
	diff_found++;
      }
      pos++;
    } while (pos < len2);
  } else if (len1 == len2 + 1) {
    do {
      if (s1[pos] != s2[pos - diff_found]) {
	if (diff_found) {
	  return 0;
	}
	diff_found++;
      }
      pos++;
    } while (pos < len1);
  } else {
    return 0;
  }
  return 1;
}

#define MAX_EQUAL_HELP_PARAMS 23

typedef struct {
  uint32_t iters_left;
  uint32_t param_ct;
  char** argv;
  uintptr_t unmatched_ct;
  uintptr_t* all_match_arr;
  uintptr_t* prefix_match_arr;
  uintptr_t* perfect_match_arr;
  uint32_t* param_lens;
  uint32_t preprint_newline;
} Help_ctrl;

void help_print(const char* cur_params, Help_ctrl* help_ctrl_ptr, uint32_t postprint_newline, const char* payload) {
  // unmatched_ct fixed during call, *unmatched_ct_ptr may decrease
  uint32_t unmatched_ct = help_ctrl_ptr->unmatched_ct;
  uint32_t print_this = 0;
  uint32_t cur_param_lens[MAX_EQUAL_HELP_PARAMS];
  char* cur_param_start[MAX_EQUAL_HELP_PARAMS];
  uint32_t arg_uidx;
  uint32_t cur_param_ct;
  uint32_t cur_param_idx;
  uint32_t arg_idx;
  uint32_t uii;
  uint32_t payload_len;
  char* payload_ptr;
  char* line_end;
  char* payload_end;
  if (help_ctrl_ptr->param_ct) {
    strcpy(g_textbuf, cur_params);
    cur_param_ct = 1;
    cur_param_start[0] = g_textbuf;
    payload_ptr = strchr(g_textbuf, '\t');
    while (payload_ptr) {
      *payload_ptr++ = '\0';
      cur_param_start[cur_param_ct++] = payload_ptr;
      payload_ptr = strchr(payload_ptr, '\t');
    }
    if (help_ctrl_ptr->iters_left) {
      if (help_ctrl_ptr->unmatched_ct) {
	arg_uidx = 0;
	if (help_ctrl_ptr->iters_left == 2) {
	  for (arg_idx = 0; arg_idx < unmatched_ct; arg_idx++) {
	    arg_uidx = next_unset_unsafe(help_ctrl_ptr->all_match_arr, arg_uidx);
	    for (cur_param_idx = 0; cur_param_idx < cur_param_ct; cur_param_idx++) {
	      if (!strcmp(cur_param_start[cur_param_idx], help_ctrl_ptr->argv[arg_uidx])) {
		SET_BIT(arg_uidx, help_ctrl_ptr->perfect_match_arr);
		SET_BIT(arg_uidx, help_ctrl_ptr->prefix_match_arr);
		SET_BIT(arg_uidx, help_ctrl_ptr->all_match_arr);
		help_ctrl_ptr->unmatched_ct -= 1;
		break;
	      }
	    }
	    arg_uidx++;
	  }
	} else {
	  for (cur_param_idx = 0; cur_param_idx < cur_param_ct; cur_param_idx++) {
	    cur_param_lens[cur_param_idx] = strlen(cur_param_start[cur_param_idx]);
	  }
	  for (arg_idx = 0; arg_idx < unmatched_ct; arg_idx++) {
	    arg_uidx = next_unset_unsafe(help_ctrl_ptr->all_match_arr, arg_uidx);
	    uii = help_ctrl_ptr->param_lens[arg_uidx];
	    for (cur_param_idx = 0; cur_param_idx < cur_param_ct; cur_param_idx++) {
	      if (cur_param_lens[cur_param_idx] > uii) {
		if (!memcmp(help_ctrl_ptr->argv[arg_uidx], cur_param_start[cur_param_idx], uii)) {
		  SET_BIT(arg_uidx, help_ctrl_ptr->prefix_match_arr);
		  SET_BIT(arg_uidx, help_ctrl_ptr->all_match_arr);
		  help_ctrl_ptr->unmatched_ct -= 1;
		  break;
		}
	      }
	    }
	    arg_uidx++;
	  }
	}
      }
    } else {
      for (cur_param_idx = 0; cur_param_idx < cur_param_ct; cur_param_idx++) {
	cur_param_lens[cur_param_idx] = strlen(cur_param_start[cur_param_idx]);
      }
      for (arg_uidx = 0; arg_uidx < help_ctrl_ptr->param_ct; arg_uidx++) {
	if (IS_SET(help_ctrl_ptr->prefix_match_arr, arg_uidx)) {
	  if (!print_this) {
	    if (IS_SET(help_ctrl_ptr->perfect_match_arr, arg_uidx)) {
	      for (cur_param_idx = 0; cur_param_idx < cur_param_ct; cur_param_idx++) {
		if (!strcmp(cur_param_start[cur_param_idx], help_ctrl_ptr->argv[arg_uidx])) {
		  print_this = 1;
		  break;
		}
	      }
	    } else {
	      uii = help_ctrl_ptr->param_lens[arg_uidx];
	      for (cur_param_idx = 0; cur_param_idx < cur_param_ct; cur_param_idx++) {
		if (cur_param_lens[cur_param_idx] > uii) {
		  if (!memcmp(help_ctrl_ptr->argv[arg_uidx], cur_param_start[cur_param_idx], uii)) {
		    print_this = 1;
		    break;
		  }
		}
	      }
	    }
	  }
	} else {
	  for (cur_param_idx = 0; cur_param_idx < cur_param_ct; cur_param_idx++) {
	    if (edit1_match(cur_param_lens[cur_param_idx], cur_param_start[cur_param_idx], help_ctrl_ptr->param_lens[arg_uidx], help_ctrl_ptr->argv[arg_uidx])) {
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
	payload_len = strlen(payload);
	if (payload[payload_len - 2] == '\n') {
	  payload_end = (char*)(&(payload[payload_len - 1]));
	} else {
	  payload_end = (char*)(&(payload[payload_len]));
	}
	if (help_ctrl_ptr->preprint_newline) {
	  putc_unlocked('\n', stdout);
	}
	help_ctrl_ptr->preprint_newline = postprint_newline;
	payload_ptr = (char*)payload;
	do {
	  line_end = strchr(payload_ptr, '\n') + 1;
	  uii = (uint32_t)(line_end - payload_ptr);
	  if (uii > 2) {
	    payload_ptr = &(payload_ptr[2]);
	    uii -= 2;
	  }
	  memcpyx(g_textbuf, payload_ptr, uii, 0);
	  fputs(g_textbuf, stdout);
	  payload_ptr = line_end;
	} while (payload_ptr < payload_end);
      }
    }
  } else {
    fputs(payload, stdout);
  }
}

int32_t disp_help(uint32_t param_ct, char** argv) {
  // yes, this is overkill.  But it should be a good template for other
  // command-line programs to use.
  uint32_t param_ctl = BITCT_TO_WORDCT(param_ct);
  int32_t retval = 0;
  Help_ctrl help_ctrl;
  uint32_t arg_uidx;
  uint32_t arg_idx;
  uint32_t net_unmatched_ct;
  int32_t col_num;
  int32_t leading_dashes;
  help_ctrl.iters_left = param_ct? 2 : 0;
  help_ctrl.param_ct = param_ct;
  help_ctrl.argv = argv;
  help_ctrl.unmatched_ct = param_ct;
  help_ctrl.param_lens = nullptr;
  help_ctrl.all_match_arr = nullptr;
  help_ctrl.argv = nullptr;
  if (param_ct) {
    help_ctrl.param_lens = (uint32_t*)malloc(param_ct * sizeof(int32_t));
    if (!help_ctrl.param_lens) {
      goto disp_help_ret_NOMEM;
    }
    help_ctrl.all_match_arr = (uintptr_t*)malloc(param_ctl * 3 * sizeof(intptr_t));
    if (!help_ctrl.all_match_arr) {
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
      help_ctrl.argv = (char**)malloc(param_ct * sizeof(char*));
      if (!help_ctrl.argv) {
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
      help_ctrl.param_lens[arg_idx] = strlen(help_ctrl.argv[arg_idx]);
    }
    fill_ulong_zero(param_ctl * 3, help_ctrl.all_match_arr);
    help_ctrl.prefix_match_arr = &(help_ctrl.all_match_arr[param_ctl]);
    help_ctrl.perfect_match_arr = &(help_ctrl.all_match_arr[param_ctl * 2]);
    help_ctrl.preprint_newline = 1;
  } else {
    help_ctrl.argv = nullptr;
    fputs(
"\nIn the command line flag definitions that follow,\n"
"  * <angle brackets> denote a required parameter, where the text between the\n"
"    angle brackets describes its nature.\n"
"  * ['square brackets + single-quotes'] denotes an optional modifier.  Use the\n"
"    EXACT text in the quotes.\n"
"  * [{bar|separated|braced|bracketed|values}] denotes a collection of mutually\n"
"    exclusive optional modifiers (again, the exact text must be used).  When\n"
"    there are no outer square brackets, one of the choices must be selected.\n"
"  * ['quoted_text='<description of value>] denotes an optional modifier that\n"
"    must begin with the quoted text, and be followed by a value with no\n"
"    whitespace in between.  '|' may also be used here to indicate mutually\n"
"    exclusive options.\n"
"  * [square brackets without quotes or braces] denote an optional parameter,\n"
"    where the text between the brackets describes its nature.\n"
"  * An ellipsis (...) indicates that you may enter multiple parameters of the\n"
"    specified type.\n"
, stdout);
    fputs(g_cmdline_format_str, stdout);
    fputs(
"Most " PROG_NAME_CAPS " runs require exactly one main input fileset.  The following flags\n"
"are available for defining its form and location:\n\n"
, stdout);
  }
  do {
    help_print("bfile\tbed\tbim\tfam", &help_ctrl, 1,
"  --bfile [prefix] : Specify .bed + .bim + .fam prefix (default '" PROG_NAME_STR "').\n"
"  --bed <filename> : Specify full name of .bed file.\n"
"  --bim <filename> : Specify full name of .bim file.\n"
"  --fam <filename> : Specify full name of .fam file.\n\n"
	       );
    help_print("file\ttfile\tlfile\tvcf\tbcf\tdata\t23file\tkeep-autoconv", &help_ctrl, 1,
"  --keep-autoconv  : With --file/--tfile/--lfile/--vcf/--bcf/--data/--23file,\n"
"                     don't delete autogenerated binary fileset at end of run.\n\n"
	       );
    help_print("file\tped\tmap", &help_ctrl, 1,
"  --file [prefix]  : Specify .ped + .map filename prefix (default '" PROG_NAME_STR "').\n"
"  --ped <filename> : Specify full name of .ped file.\n"
"  --map <filename> : Specify full name of .map file.\n\n"
	       );
    help_print("bfile\tfam\tfile\tped\tno-fid\tno-parents\tno-sex\tno-pheno", &help_ctrl, 1,
"  --no-fid         : .fam/.ped file does not contain column 1 (family ID).\n"
"  --no-parents     : .fam/.ped file does not contain columns 3-4 (parents).\n"
"  --no-sex         : .fam/.ped file does not contain column 5 (sex).\n"
"  --no-pheno       : .fam/.ped file does not contain column 6 (phenotype).\n\n"
	       );
    help_print("tfile\ttped\ttfam", &help_ctrl, 1,
"  --tfile [prefix] : Specify .tped + .tfam filename prefix (default '" PROG_NAME_STR "').\n"
"  --tped <fname>   : Specify full name of .tped file.\n"
"  --tfam <fname>   : Specify full name of .tfam file.\n\n"
	       );
    help_print("lfile\treference\tallele-count", &help_ctrl, 1,
"  --lfile [prefix] : Specify .lgen + .map + .fam (long-format fileset) prefix.\n"
"  --lgen <fname>   : Specify full name of .lgen file.\n"
"  --reference <fn> : Specify default allele file accompanying .lgen input.\n"
"  --allele-count   : When used with --lfile/--lgen + --reference, specifies\n"
"                     that the .lgen file contains reference allele counts.\n\n"
	       );
    help_print("vcf\tbcf", &help_ctrl, 1,
"  --vcf <filename> : Specify full name of .vcf or .vcf.gz file.\n"
"  --bcf <filename> : Specify full name of BCF2 file.\n\n"
	       );
    help_print("data\tgen\tbgen\tsample", &help_ctrl, 1,
"  --data [prefix]  : Specify Oxford .gen + .sample prefix (default '" PROG_NAME_STR "').\n"
"  --gen <filename> : Specify full name of .gen or .gen.gz file.\n"
"  --bgen <f> ['snpid-chr'] : Specify full name of .bgen file.\n"
"  --sample <fname> : Specify full name of .sample file.\n\n"
    	       );
    help_print("23file", &help_ctrl, 1,
"  --23file <fname> [FID] [IID] [sex] [pheno] [pat. ID] [mat. ID] :\n"
"    Specify 23andMe input file.\n\n"
	       );
#ifndef STABLE_BUILD
    help_print("cfile\tcnv-list\tgfile", &help_ctrl, 1,
"  --cfile <prefix> : Specify .cnv + .fam + .cnv.map (segmental CNV) prefix.\n"
"  --cnv-list <fn>  : Specify full name of .cnv file.\n"
"  --gfile <prefix> : Specify .gvar + .fam + .map (genetic variant) prefix.\n\n"
	       );
#endif
    help_print("grm\tgrm-gz\tgrm-bin\trel-cutoff\tgrm-cutoff", &help_ctrl, 1,
"  --grm-gz [prfx]  : Specify .grm.gz + .grm.id (GCTA rel. matrix) prefix.\n"
"  --grm-bin [prfx] : Specify .grm.bin + .grm.N.bin + .grm.id (GCTA triangular\n"
"                     binary relationship matrix) filename prefix.\n\n"
	       );
    help_print("dummy", &help_ctrl, 1,
"  --dummy <sample ct> <SNP ct> [missing geno freq] [missing pheno freq]\n"
"          [{acgt | 1234 | 12}] ['scalar-pheno']\n"
"    This generates a fake input dataset with the specified number of samples\n"
"    and SNPs.  By default, the missing genotype and phenotype frequencies are\n"
"    zero, and genotypes are As and Bs (change the latter with\n"
"    'acgt'/'1234'/'12').  The 'scalar-pheno' modifier causes a normally\n"
"    distributed scalar phenotype to be generated instead of a binary one.\n\n"
	       );
    help_print("simulate\tsimulate-qt", &help_ctrl, 1,
"  --simulate <simulation parameter file> [{tags | haps}] [{acgt | 1234 | 12}]\n"
"  --simulate-qt <sim. parameter file> [{tags | haps}] [{acgt | 1234 | 12}]\n"
"    --simulate generates a fake input dataset with disease-associated SNPs,\n"
"    while --simulate-qt generates a dataset with quantitative trait loci.\n\n"
	       );
    if (!param_ct) {
      fputs(
"Output files have names of the form '" PROG_NAME_STR ".<extension>' by default.  You can\n"
"change the '" PROG_NAME_STR "' prefix with\n\n"
, stdout);
    }
    help_print("out", &help_ctrl, 1,
"  --out <prefix>   : Specify prefix for output files.\n\n"
	       );
    if (!param_ct) {
      fputs(
"Most runs also require at least one of the following commands:\n\n"
, stdout);
    }
    help_print("make-bed", &help_ctrl, 1,
"  --make-bed\n"
"    Create a new binary fileset.  Unlike the automatic text-to-binary\n"
"    converters (which only heed chromosome filters), this supports all of\n"
"    " PROG_NAME_CAPS "'s filtering flags.\n"
	       );
    help_print("make-just-bim\tmake-just-fam", &help_ctrl, 1,
"  --make-just-bim\n"
"  --make-just-fam\n"
"    Variants of --make-bed which only write a new .bim or .fam file.  Can be\n"
"    used with only .bim/.fam input.\n"
"    USE THESE CAUTIOUSLY.  It is very easy to desynchronize your binary\n"
"    genotype data and your .bim/.fam indexes if you use these commands\n"
"    improperly.  If you have any doubt, stick with --make-bed.\n\n"
	       );

    help_print("recode\trecode12\ttab\ttranspose\trecode-lgen\trecodeAD\trecodead\trecodeA\trecodea\trecode-rlist\trecode-allele\tlist\twith-reference\trecode-vcf\tfid\tiid\trecode-beagle\trecode-bimbam\trecode-fastphase\trecodeHV\trecodehv\trecode-structure\texport", &help_ctrl, 1,
"  --recode <output format> [{01 | 12}] [{tab | tabx | spacex | bgz | gen-gz}]\n"
"           ['include-alt'] ['omit-nonmale-y']\n"
"    Create a new text fileset with all filters applied.  The following output\n"
"    formats are supported:\n"
"    * '23': 23andMe 4-column format.  This can only be used on a single\n"
"      sample's data (--keep may be handy), and does not support multicharacter\n"
"      allele codes.\n"
"    * 'A': Sample-major additive (0/1/2) coding, suitable for loading from R.\n"
"      If you need uncounted alleles to be named in the header line, add the\n"
"      'include-alt' modifier.\n"
"    * 'AD': Sample-major additive (0/1/2) + dominant (het=1/hom=0) coding.\n"
"      Also supports 'include-alt'.\n"
"    * 'A-transpose': Variant-major 0/1/2.\n"
"    * 'beagle': Unphased per-autosome .dat and .map files, readable by early\n"
"      BEAGLE versions.\n"
"    * 'beagle-nomap': Single .beagle.dat file.\n"
"    * 'bimbam': Regular BIMBAM format.\n"
"    * 'bimbam-1chr': BIMBAM format, with a two-column .pos.txt file.  Does not\n"
"      support multiple chromosomes.\n"
"    * 'fastphase': Per-chromosome fastPHASE files, with\n"
"      .chr-<chr #>.recode.phase.inp filename extensions.\n"
"    * 'fastphase-1chr': Single .recode.phase.inp file.  Does not support\n"
"      multiple chromosomes.\n"
"    * 'HV': Per-chromosome Haploview files, with .chr-<chr #>{.ped,.info}\n"
"      filename extensions.\n"
"    * 'HV-1chr': Single Haploview .ped + .info file pair.  Does not support\n"
"      multiple chromosomes.\n"
"    * 'lgen': PLINK 1 long-format (.lgen + .fam + .map), loadable with --lfile.\n"
"    * 'lgen-ref': .lgen + .fam + .map + .ref, loadable with --lfile +\n"
"       --reference.\n"
"    * 'list': Single genotype-based list, up to 4 lines per variant.  To omit\n"
"      nonmale genotypes on the Y chromosome, add the 'omit-nonmale-y' modifier.\n"
"    * 'rlist': .rlist + .fam + .map fileset, where the .rlist file is a\n"
"      genotype-based list which omits the most common genotype for each\n"
"      variant.  Also supports 'omit-nonmale-y'.\n"
"    * 'oxford': Oxford-format .gen + .sample.  With the 'gen-gz' modifier, the\n"
"      .gen file is gzipped.\n"
"    * 'ped': PLINK 1 sample-major (.ped + .map), loadable with --file.\n"
"    * 'compound-genotypes': Same as 'ped', except that the space between each\n"
"      pair of same-variant allele codes is removed.\n"
"    * 'structure': Structure-format.\n"
"    * 'transpose': PLINK 1 variant-major (.tped + .tfam), loadable with\n"
"      --tfile.\n"
"    * 'vcf', 'vcf-fid', 'vcf-iid': VCFv4.2.  'vcf-fid' and 'vcf-iid' cause\n"
"      family IDs or within-family IDs respectively to be used for the sample\n"
"      IDs in the last header row, while 'vcf' merges both IDs and puts an\n"
"      underscore between them.  If the 'bgz' modifier is added, the VCF file is\n"
"      block-gzipped.\n"
"      The A2 allele is saved as the reference and normally flagged as not based\n"
"      on a real reference genome (INFO:PR).  When it is important for reference\n"
"      alleles to be correct, you'll also want to include --a2-allele and\n"
"      --real-ref-alleles in your command.\n"
"    In addition,\n"
"    * The '12' modifier causes A1 (usually minor) alleles to be coded as '1'\n"
"      and A2 alleles to be coded as '2', while '01' maps A1 -> 0 and A2 -> 1.\n"
"    * The 'tab' modifier makes the output mostly tab-delimited instead of\n"
"      mostly space-delimited.  'tabx' and 'spacex' force all tabs and all\n"
"      spaces, respectively.\n\n"
	       );
    help_print("flip-scan\tflip-scan-verbose\tflipscan", &help_ctrl, 1,
"  --flip-scan ['verbose']\n"
"    (alias: --flipscan)\n"
"    LD-based scan for case/control strand inconsistency.\n\n"
	       );
    help_print("write-covar", &help_ctrl, 1,
"  --write-covar\n"
"    If a --covar file is loaded, --make-bed/--make-just-fam and --recode\n"
"    automatically generate an updated version (with all filters applied).\n"
"    However, if you do not wish to simultaneously generate a new genotype file,\n"
"    you can use --write-covar to just produce a pruned covariate file.\n\n"
	       );
    help_print("write-cluster", &help_ctrl, 1,
"  --write-cluster ['omit-unassigned']\n"
"    If clusters are specified with --within/--family, this generates a new\n"
"    cluster file (with all filters applied).  The 'omit-unassigned' modifier\n"
"    causes unclustered samples to be omitted from the file; otherwise their\n"
"    cluster is 'NA'.\n\n"
	       );
    help_print("write-set\tset-table", &help_ctrl, 1,
"  --write-set\n"
"  --set-table\n"
"    If sets have been defined, --write-set dumps 'END'-terminated set\n"
"    membership lists to <output prefix>.set, while --set-table writes a\n"
"    variant-by-set membership table to <output prefix>.set.table.\n\n"
	       );
    help_print("merge\tbmerge\tmerge-list\tmerge-mode", &help_ctrl, 1,
"  --merge <.ped filename> <.map filename>\n"
"  --merge <text fileset prefix>\n"
"  --bmerge <.bed filename> <.bim filename> <.fam filename>\n"
"  --bmerge <binary fileset prefix>\n"
"    Merge the given fileset with the initially loaded fileset, writing the\n"
"    result to <output prefix>.bed + .bim + .fam.  (It is no longer necessary to\n"
"    simultaneously specify --make-bed.)\n"
"  --merge-list <filename>\n"
"    Merge all filesets named in the text file with the reference fileset, if\n"
"    one was specified.  (However, this can also be used *without* a reference;\n"
"    in that case, the newly created fileset is then treated as the reference by\n"
"    most other PLINK operations.)  The text file is interpreted as follows:\n"
"    * If a line contains only one name, it is assumed to be the prefix for a\n"
"      binary fileset.\n"
"    * If a line contains exactly two names, they are assumed to be the full\n"
"      filenames for a text fileset (.ped first, then .map).\n"
"    * If a line contains exactly three names, they are assumed to be the full\n"
"      filenames for a binary fileset (.bed, then .bim, then .fam).\n\n"
	       );
    help_print("write-snplist\tlist-23-indels", &help_ctrl, 1,
"  --write-snplist\n"
"  --list-23-indels\n"
"    --write-snplist writes a .snplist file listing the names of all variants\n"
"    which pass the filters and inclusion thresholds you've specified, while\n"
"    --list-23-indels writes the subset with 23andMe-style indel calls (D/I\n"
"    allele codes).\n\n"
	       );
    help_print("list-duplicate-vars", &help_ctrl, 1,
"  --list-duplicate-vars ['require-same-ref'] ['ids-only'] ['suppress-first']\n"
"    --list-duplicate-vars writes a .dupvar file describing all groups of\n"
"    variants with matching positions and allele codes.\n"
"    * By default, A1/A2 allele assignments are ignored; use 'require-same-ref'\n"
"      to override this.\n"
"    * Normally, the report contains position and allele codes.  To remove them\n"
"      (and produce a file directly usable with e.g. --extract/--exclude), use\n"
"      'ids-only'.  Note that this command will fail in 'ids-only' mode if any\n"
"      of the reported IDs are not unique.\n"
"    * 'suppress-first' causes the first variant ID in each group to be omitted\n"
"      from the report.\n\n"
	       );
    help_print("freq\tfreqx\tfrqx\tcounts", &help_ctrl, 1,
"  --freq [{counts | case-control}] ['gz']\n"
"  --freqx ['gz']\n"
"    --freq generates a basic allele frequency (or count, if the 'counts'\n"
"    modifier is present) report.  This can be combined with --within/--family\n"
"    to produce a cluster-stratified allele frequency/count report instead, or\n"
"    the 'case-control' modifier to report case and control allele frequencies\n"
"    separately.\n"
"    --freqx generates a more detailed genotype count report, designed for use\n"
"    with --read-freq.\n\n"
		);
    help_print("missing", &help_ctrl, 1,
"  --missing ['gz']\n"
"    Generate sample- and variant-based missing data reports.  If clusters are\n"
"    defined, the variant-based report is cluster-stratified.  'gz' causes the\n"
"    output files to be gzipped.\n"
"    Unlike most other commands, this doesn't treat het. haploids as missing.\n\n"
	       );
    help_print("test-mishap", &help_ctrl, 1,
"  --test-mishap\n"
"    Check for association between missing calls and flanking haplotypes.\n\n"
               );
    help_print("hardy\thardy2", &help_ctrl, 1,
"  --hardy ['midp'] ['gz']\n"
"    Generate a Hardy-Weinberg exact test p-value report.  (This does NOT\n"
"    simultaneously filter on the p-value any more; use --hwe for that.)  With\n"
"    the 'midp' modifier, the test applies the mid-p adjustment described in\n"
"    Graffelman J, Moreno V (2013) The mid p-value in exact tests for\n"
"    Hardy-Weinberg Equilibrium.\n\n"
	       );
    help_print("mendel", &help_ctrl, 1,
"  --mendel ['summaries-only']\n"
"    Generate a Mendel error report.  The 'summaries-only' modifier causes the\n"
"    .mendel file (listing every single error) to be skipped.\n\n"
	       );
    help_print("het\tibc", &help_ctrl, 1,
"  --het ['small-sample'] ['gz']\n"
"  --ibc\n"
"    Estimate inbreeding coefficients.  --het reports method-of-moments\n"
"    estimates, while --ibc calculates all three values described in Yang J, Lee\n"
"    SH, Goddard ME and Visscher PM (2011) GCTA: A Tool for Genome-wide Complex\n"
"    Trait Analysis.  (That paper also describes the relationship matrix\n"
"    computation we reimplement.)\n"
"    * These functions require decent MAF estimates.  If there are very few\n"
"      samples in your immediate fileset, --read-freq is practically mandatory\n"
"      since imputed MAFs are wildly inaccurate in that case.\n"
"    * They also assume the marker set is in approximate linkage equilibrium.\n"
"    * By default, --het omits the n/(n-1) multiplier in Nei's expected\n"
"      homozygosity formula.  The 'small-sample' modifier causes it to be\n"
"      included, while forcing --het to use MAFs imputed from founders in the\n"
"      immediate dataset.\n\n"
	       );
    help_print("check-sex\timpute-sex\tupdate-sex\tsex-check", &help_ctrl, 1,
"  --check-sex [female max F] [male min F]\n"
"  --check-sex ycount [female max F] [male min F] [female max Y obs]\n"
"                     [male min Y obs]\n"
"  --check-sex y-only [female max Y obs] [male min Y obs]\n"
"  --impute-sex [female max F] [male min F]\n"
"  --impute-sex ycount [female max F] [male min F] [female max Y obs]\n"
"                      [male min Y obs]\n"
"  --impute-sex y-only [female max Y obs] [male min Y obs]\n"
"    --check-sex normally compares sex assignments in the input dataset with\n"
"    those imputed from X chromosome inbreeding coefficients.\n"
"    * Make sure that the X chromosome pseudo-autosomal region has been split\n"
"      off (with e.g. --split-x) before using this.\n"
"    * You also need decent MAF estimates (so, with very few samples in your\n"
"      immediate fileset, use --read-freq), and your marker set should be in\n"
"      approximate linkage equilibrium.\n"
"    * By default, F estimates smaller than 0.2 yield female calls, and values\n"
"      larger than 0.8 yield male calls.  If you pass numeric parameter(s) to\n"
"      --check-sex, the first two control these thresholds.\n"
"    There are now two modes which consider Y chromosome data.\n"
"    * In 'ycount' mode, gender is still imputed from the X chromosome, but\n"
"      female calls are downgraded to ambiguous whenever more than 0 nonmissing\n"
"      Y genotypes are present, and male calls are downgraded when fewer than 0\n"
"      are present.  (Note that these are counts, not rates.)  These thresholds\n"
"      are controllable with --check-sex ycount's optional 3rd and 4th numeric\n"
"      parameters.\n"
"    * In 'y-only' mode, gender is imputed from nonmissing Y genotype counts.\n"
"      The male minimum threshold defaults to 1 instead of zero in this case.\n"
"    --impute-sex changes sex assignments to the imputed values, and is\n"
"    otherwise identical to --check-sex.  It must be used with\n"
"    --make-bed/--recode/--write-covar.\n\n"
	       );
    help_print("fst\tFst", &help_ctrl, 1,
"  --fst ['case-control']\n"
"    (alias: --Fst)\n"
"    Estimate Wright's Fst for each autosomal diploid variant using the method\n"
"    introduced in Weir BS, Cockerham CC (1984) Estimating F-statistics for the\n"
"    analysis of population structure, given a set of subpopulations defined via\n"
"    --within.  Raw and weighted global means are also reported.\n"
"    * If you're interested in the global means, it is usually best to perform\n"
"      this calculation on a marker set in approximate linkage equilibrium.\n"
"    * If you have only two subpopulations, you can represent them with\n"
"      case/control status and use the 'case-control' modifier.\n\n"
	       );
    help_print("indep\tindep-pairwise\tindep-pairphase", &help_ctrl, 1,
"  --indep <window size>['kb'] <step size (variant ct)> <VIF threshold>\n"
"  --indep-pairwise <window size>['kb'] <step size (variant ct)> <r^2 threshold>\n"
"  --indep-pairphase <window size>['kb'] <step size (variant ct)> <r^2 thresh>\n"
"    Generate a list of markers in approximate linkage equilibrium.  With the\n"
"    'kb' modifier, the window size is in kilobase instead of variant count\n"
"    units.  (Pre-'kb' space is optional, i.e. \"--indep-pairwise 500 kb 5 0.5\"\n"
"    and \"--indep-pairwise 500kb 5 0.5\" have the same effect.)\n"
"    Note that you need to rerun " PROG_NAME_CAPS " using --extract or --exclude on the\n"
"    .prune.in/.prune.out file to apply the list to another computation.\n\n"
		);
    help_print("r\tr2\tmatrix\tinter-chr\tD\tdprime\twith-freqs\tld", &help_ctrl, 1,
"  --r [{square | square0 | triangle | inter-chr}] [{gz | bin | bin4}]\n"
"      ['spaces'] ['in-phase'] [{d | dprime | dprime-signed}] ['with-freqs']\n"
"      ['yes-really']\n"
"  --r2 [{square | square0 | triangle | inter-chr}] [{gz | bin | bin4}]\n"
"       ['spaces'] ['in-phase'] [{d | dprime | dprime-signed}] ['with-freqs']\n"
"       ['yes-really']\n"
"    LD statistic reports.  --r yields raw inter-variant correlations, while\n"
"    --r2 reports their squares.  You can request results for all pairs in\n"
"    matrix format (if you specify 'bin' or one of the shape modifiers), all\n"
"    pairs in table format ('inter-chr'), or a limited window in table format\n"
"    (default).\n"
"    * The 'gz' modifier causes the output text file to be gzipped.\n"
"    * 'bin' causes the output matrix to be written in double-precision binary\n"
"      format, while 'bin4' specifics single-precision binary.  The matrix is\n"
"      square if no shape is explicitly specified.\n"
"    * By default, text matrices are tab-delimited; 'spaces' switches this.\n"
"    * 'in-phase' adds a column with in-phase allele pairs to table-formatted\n"
"      reports.  (This cannot be used with very long allele codes.)\n"
"    * 'dprime' adds the absolute value of Lewontin's D-prime statistic to\n"
"      table-formatted reports, and forces both r/r^2 and D-prime to be based on\n"
"      the maximum likelihood solution to the cubic equation discussed in Gaunt\n"
"      T, Rodriguez S, Day I (2007) Cubic exact solutions for the estimation of\n"
"      pairwise haplotype frequencies.\n"
"      'dprime-signed' keeps the sign, while 'd' skips division by D_{max}.\n"
"    * 'with-freqs' adds MAF columns to table-formatted reports.\n"
"    * Since the resulting file can easily be huge, you're required to add the\n"
"      'yes-really' modifier when requesting an unfiltered, non-distributed all\n"
"      pairs computation on more than 400k variants.\n"
"    * These computations can be subdivided with --parallel (even when the\n"
"      'square' modifier is active).\n"
"  --ld <variant ID> <variant ID> ['hwe-midp']\n"
"    This displays haplotype frequencies, r^2, and D' for a single pair of\n"
"    variants.  When there are multiple biologically possible solutions to the\n"
"    haplotype frequency cubic equation, all are displayed (instead of just the\n"
"    maximum likelihood solution identified by --r/--r2), along with HWE exact\n"
"    test statistics.\n\n"
	       );
    help_print("show-tags", &help_ctrl, 1,
"  --show-tags <filename>\n"
"  --show-tags all\n"
"    * If a file is specified, list all variants which tag at least one variant\n"
"      named in the file.  (This will normally be a superset of the original\n"
"      list, since a variant is considered to tag itself here.)\n"
"    * If 'all' mode is specified, for each variant, each *other* variant which\n"
"      tags it is reported.\n\n"
	       );
    help_print("blocks\thap\thap-all\thap-assoc\thap-freq\thap-impute\thap-impute-verbose\thap-linear\thap-logistic\thap-max-phase\thap-min-phase-prob\thap-miss\thap-omnibus\thap-only\thap-phase\thap-phase-wide\thap-pp\thap-snps\thap-tdt\thap-window\tchap\twhap", &help_ctrl, 1,
"  --blocks ['no-pheno-req'] ['no-small-max-span']\n"
"    Estimate haplotype blocks, via Haploview's interpretation of the block\n"
"    definition suggested by Gabriel S et al. (2002) The Structure of Haplotype\n"
"    Blocks in the Human Genome.\n"
"    * Normally, samples with missing phenotypes are not considered by this\n"
"      computation; the 'no-pheno-req' modifier lifts this restriction.\n"
"    * Normally, size-2 blocks may not span more than 20kb, and size-3 blocks\n"
"      are limited to 30kb.  The 'no-small-max-span' modifier removes these\n"
"      limits.\n"
"    The .blocks file is valid input for PLINK 1.07's --hap command.  However,\n"
"    the --hap... family of flags has not been reimplemented in PLINK 1.9 due to\n"
"    poor phasing accuracy relative to other software; for now, we recommend\n"
"    using BEAGLE instead of PLINK for case/control haplotype association\n"
"    analysis.  (You can use \"--recode beagle\" to export data to BEAGLE 3.3.)\n"
"    We apologize for the inconvenience, and plan to develop variants of the\n"
"    --hap... flags which handle pre-phased data effectively.\n\n"
	       );
    help_print("distance", &help_ctrl, 1,
"  --distance [{square | square0 | triangle}] [{gz | bin | bin4}] ['ibs']\n"
"             ['1-ibs'] ['allele-ct'] ['flat-missing']\n"
"    Write a lower-triangular tab-delimited table of (weighted) genomic\n"
"    distances in allele count units to <output prefix>.dist, and a list of the\n"
"    corresponding sample IDs to <output prefix>.dist.id.  The first row of the\n"
"    .dist file contains a single <genome 1-genome 2> distance, the second row\n"
"    has the <genome 1-genome 3> and <genome 2-genome 3> distances in that\n"
"    order, etc.\n"
"    * It is usually best to perform this calculation on a marker set in\n"
"      approximate linkage equilibrium.\n"
"    * If the 'square' or 'square0' modifier is present, a square matrix is\n"
"      written instead; 'square0' fills the upper right triangle with zeroes.\n"
"    * If the 'gz' modifier is present, a compressed .dist.gz file is written\n"
"      instead of a plain text file.\n"
"    * If the 'bin' modifier is present, a binary (square) matrix of\n"
"      double-precision floating point values, suitable for loading from R, is\n"
"      instead written to <output prefix>.dist.bin.  ('bin4' specifies\n"
"      single-precision numbers instead.)  This can be combined with 'square0'\n"
"      if you still want the upper right zeroed out, or 'triangle' if you don't\n"
"      want to pad the upper right at all.\n"
"    * If the 'ibs' modifier is present, an identity-by-state matrix is written\n"
"      to <output prefix>.mibs.  '1-ibs' causes distances expressed as genomic\n"
"      proportions (i.e. 1 - IBS) to be written to <output prefix>.mdist.\n"
"      Combine with 'allele-ct' if you want to generate the usual .dist file as\n"
"      well.\n"
"    * By default, distance rescaling in the presence of missing genotype calls\n"
"      is sensitive to allele count distributions: if variant A contributes, on\n"
"      average, twice as much to other pairwise distances as variant B, a\n"
"      missing call at variant A will result in twice as large of a missingness\n"
"      correction.  To turn this off (because e.g. your missing calls are highly\n"
"      nonrandom), use the 'flat-missing' modifier.\n"
"    * The computation can be subdivided with --parallel.\n"
	       );
    help_print("distance-matrix\tibs-matrix\tmatrix", &help_ctrl, 1,
"  --distance-matrix\n"
"  --ibs-matrix\n"
"    These deprecated commands are equivalent to \"--distance 1-ibs flat-missing\n"
"    square\" and \"--distance ibs flat-missing square\", respectively, except that\n"
"    they generate space- instead of tab-delimited text matrices.\n\n"
		);
    help_print("make-rel", &help_ctrl, 1,
"  --make-rel [{square | square0 | triangle}] [{gz | bin | bin4}]\n"
"             [{cov | ibc2 | ibc3}]\n"
"    Write a lower-triangular variance-standardized realized relationship matrix\n"
"    to <output prefix>.rel, and corresponding IDs to <output prefix>.rel.id.\n"
"    * It is usually best to perform this calculation on a marker set in\n"
"      approximate linkage equilibrium.\n"
"    * 'square', 'square0', 'triangle', 'gz', 'bin', and 'bin4' act as they do\n"
"      on --distance.\n"
"    * The 'cov' modifier removes the variance standardization step, causing a\n"
"      covariance matrix to be calculated instead.\n"
"    * By default, the diagonal elements in the relationship matrix are based on\n"
"      --ibc's Fhat1; use the 'ibc2' or 'ibc3' modifiers to base them on Fhat2\n"
"      or Fhat3 instead.\n"
"    * The computation can be subdivided with --parallel.\n"
               );
    help_print("make-grm\tmake-grm-bin\tgrm\tgrm-bin\tmake-grm-gz", &help_ctrl, 1,
"  --make-grm-gz ['no-gz'] [{cov | ibc2 | ibc3}]\n"
"  --make-grm-bin [{cov | ibc2 | ibc3}]\n"
"    --make-grm-gz writes the relationships in GCTA's original gzipped list\n"
"    format, which describes one pair per line, while --make-grm-bin writes them\n"
"    in GCTA 1.1+'s single-precision triangular binary format.  Note that these\n"
"    formats explicitly report the number of valid observations (where neither\n"
"    sample has a missing call) for each pair, which is useful input for some\n"
"    scripts.\n"
"    These computations can be subdivided with --parallel.\n\n"
	       );
    help_print("rel-cutoff\tgrm-cutoff", &help_ctrl, 1,
"  --rel-cutoff [val]\n"
"    (alias: --grm-cutoff)\n"
"    Exclude one member of each pair of samples with relatedness greater than\n"
"    the given cutoff value (default 0.025).  If no later operation will cause\n"
"    the list of remaining samples to be written to disk, this will save it to\n"
"    <output prefix>.rel.id.\n"
"    Note that maximizing the remaining sample size is equivalent to the NP-hard\n"
"    maximum independent set problem, so we use a greedy algorithm instead of\n"
"    guaranteeing optimality.  (Use the --make-rel and --keep/--remove flags if\n"
"    you want to try to do better.)\n\n"
	       );
    help_print("ibs-test\tgroupdist", &help_ctrl, 1,
"  --ibs-test [permutation count]\n"
"  --groupdist [iters] [d]\n"
"    Given case/control phenotype data, these commands consider three subsets of\n"
"    the distance matrix: pairs of affected samples, affected-unaffected pairs,\n"
"    and pairs of unaffected samples.  Each of these subsets has a distribution\n"
"    of pairwise genomic distances; --ibs-test uses permutation to estimate\n"
"    p-values re: which types of pairs are most similar, while --groupdist\n"
"    focuses on the differences between the centers of these distributions and\n"
"    estimates standard errors via delete-d jackknife.\n\n"
	       );
    help_print("regress-distance\tregress-rel", &help_ctrl, 1,
"  --regress-distance [iters] [d]\n"
"    Linear regression of pairwise genomic distances on pairwise average\n"
"    phenotypes and vice versa, using delete-d jackknife for standard errors.  A\n"
"    scalar phenotype is required.\n"
"    * With less than two parameters, d is set to <number of people>^0.6 rounded\n"
"      down.  With no parameters, 100k iterations are run.\n"
"  --regress-rel [iters] [d]\n"
"    Linear regression of pairwise genomic relationships on pairwise average\n"
"    phenotypes, and vice versa.  Defaults for iters and d are the same as for\n"
"    --regress-distance.\n\n"
	       );
    help_print("genome\tZ-genome\trel-check\timpossible\tnudge\tgenome-full\tunbounded", &help_ctrl, 1,
"  --genome ['gz'] ['rel-check'] ['full'] ['unbounded'] ['nudge']\n"
"    Generate an identity-by-descent report.\n"
"    * It is usually best to perform this calculation on a marker set in\n"
"      approximate linkage equilibrium.\n"
"    * The 'rel-check' modifier excludes pairs of samples with different FIDs\n"
"      from the final report.\n"
"    * 'full' adds raw pairwise comparison data to the report.\n"
"    * The P(IBD=0/1/2) estimator employed by this command sometimes yields\n"
"      numbers outside the range [0,1]; by default, these are clipped.  The\n"
"      'unbounded' modifier turns off this clipping.\n"
"    * Then, when PI_HAT^2 < P(IBD=2), 'nudge' adjusts the final P(IBD=0/1/2)\n"
"      estimates to a theoretically possible configuration.\n"
"    * The computation can be subdivided with --parallel.\n\n"
		);
    help_print("homozyg\thomozyg-snp\thomozyg-kb\thomozyg-density\thomozyg-gap\thomozyg-het\thomozyg-window-snp\thomozyg-window-het\thomozyg-window-missing\thomozyg-window-threshold", &help_ctrl, 1,
"  --homozyg [{group | group-verbose}] ['consensus-match'] ['extend']\n"
"            ['subtract-1-from-lengths']\n"
"  --homozyg-snp <min var count>\n"
"  --homozyg-kb <min length>\n"
"  --homozyg-density <max inverse density (kb/var)>\n"
"  --homozyg-gap <max internal gap kb length>\n"
"  --homozyg-het <max hets>\n"
"  --homozyg-window-snp <scanning window size>\n"
"  --homozyg-window-het <max hets in scanning window hit>\n"
"  --homozyg-window-missing <max missing calls in scanning window hit>\n"
"  --homozyg-window-threshold <min scanning window hit rate>\n"
"    These commands request a set of run-of-homozygosity reports, and allow you\n"
"    to customize how they are generated.\n"
"    * If you're satisfied with all the default settings described below, just\n"
"      use --homozyg with no modifiers.  Otherwise, --homozyg lets you change a\n"
"      few binary settings:\n"
"      * 'group[-verbose]' adds a report on pools of overlapping runs of\n"
"        homozygosity.  (Automatically set when --homozyg-match is present.)\n"
"      * With 'group[-verbose]', 'consensus-match' causes pairwise segmental\n"
"        matches to be called based on the variants in the pool's consensus\n"
"        segment, rather than the variants in the pairwise intersection.\n"
"      * Due to how the scanning window algorithm works, it is possible for a\n"
"        reported ROH to be adjacent to a few homozygous variants.  The 'extend'\n"
"        modifier causes them to be included in the reported ROH if that\n"
"        wouldn't cause a violation of the --homozyg-density bound.\n"
"      * By default, segment bp lengths are calculated as <end bp position> -\n"
"        <start bp position> + 1.  Therefore, reports normally differ slightly\n"
"        from PLINK 1.07, which does not add 1 at the end.  For testing\n"
"        purposes, you can use the 'subtract-1-from-lengths' modifier to apply\n"
"        the old formula.\n"
"    * By default, only runs of homozygosity containing at least 100 variants,\n"
"      and of total length >= 1000 kilobases, are noted.  You can change these\n"
"      minimums with --homozyg-snp and --homozyg-kb, respectively.\n"
"    * By default, a ROH must have at least one variant per 50 kb on average;\n"
"      change this bound with --homozyg-density.\n"
"    * By default, if two consecutive variants are more than 1000 kb apart, they\n"
"      cannot be in the same ROH; change this bound with --homozyg-gap.\n"
"    * By default, a ROH can contain an unlimited number of heterozygous calls;\n"
"      you can impose a limit with --homozyg-het.\n"
"    * By default, the scanning window contains 50 variants; change this with\n"
"      --homozyg-window-snp.\n"
"    * By default, a scanning window hit can contain at most 1 heterozygous\n"
"      call and 5 missing calls; change these limits with --homozyg-window-het\n"
"      and --homozyg-window-missing, respectively.\n"
"    * By default, for a variant to be eligible for inclusion in a ROH, the hit\n"
"      rate of all scanning windows containing the variant must be at least\n"
"      0.05; change this threshold with --homozyg-window-threshold.\n\n"
	       );
    help_print("cluster\tcc\tgroup-avg\tgroup-average\tcluster-missing", &help_ctrl, 1,
"  --cluster ['cc'] [{group-avg | old-tiebreaks}] ['missing'] ['only2']\n"
"    Cluster samples using a pairwise similarity statistic (normally IBS).\n"
"    * The 'cc' modifier forces every cluster to have at least one case and one\n"
"      control.\n"
"    * The 'group-avg' modifier causes clusters to be joined based on average\n"
"      instead of minimum pairwise similarity.\n"
"    * The 'missing' modifier causes clustering to be based on\n"
"      identity-by-missingness instead of identity-by-state, and writes a\n"
"      space-delimited identity-by-missingness matrix to disk.\n"
"    * The 'only2' modifier causes only a .cluster2 file (which is valid input\n"
"      for --within) to be written; otherwise 2 other files will be produced.\n"
"    * By default, IBS ties are not broken in the same manner as PLINK 1.07, so\n"
"      final cluster solutions tend to differ.  This is generally harmless.\n"
"      However, to simplify testing, you can use the 'old-tiebreaks' modifier to\n"
"      force emulation of the old algorithm.\n\n"
	       );
#ifndef NOLAPACK
    help_print("pca\tmake-rel\tmake-grm\tmake-grm-gz\tmake-grm-bin", &help_ctrl, 1,
"  --pca [count] ['header'] ['tabs'] ['var-wts']\n"
"    Calculates a variance-standardized relationship matrix (use\n"
"    --make-rel/--make-grm-gz/--make-grm-bin to dump it), and extracts the top\n"
"    20 principal components.\n"
"    * It is usually best to perform this calculation on a marker set in\n"
"      approximate linkage equilibrium.\n"
"    * You can change the number of PCs by passing a numeric parameter.\n"
"    * The 'header' modifier adds a header line to the .eigenvec output file.\n"
"      (For compatibility with the GCTA flag of the same name, the default is no\n"
"      header line.)\n"
"    * The 'tabs' modifier causes the .eigenvec file(s) to be tab-delimited.\n"
"    * The 'var-wts' modifier requests an additional .eigenvec.var file with PCs\n"
"      expressed as variant weights instead of sample weights.\n\n"
	       );
#endif
    help_print("neighbour\tneighbor", &help_ctrl, 1,
"  --neighbour <n1> <n2>\n"
"    (alias: --neighbor)\n"
"    Report IBS distances from each sample to their n1th- to n2th-nearest\n"
"    neighbors, associated Z-scores, and the identities of those neighbors.\n"
"    Useful for outlier detection.\n\n"
	       );
    help_print("assoc\tmodel\tfisher\tperm\tmperm\tperm-count\tcounts\tp2\tset-test\tmodel-dom\tmodel-gen\tmodel-rec\tmodel-trend\tgenedrop\tqt-means\ttrend", &help_ctrl, 1,
"  --assoc ['perm' | 'mperm='<value>] ['perm-count'] [{fisher | fisher-midp}]\n"
"          ['counts'] ['set-test']\n"
"  --assoc ['perm' | 'mperm='<value>] ['perm-count'] ['qt-means'] ['lin']\n"
"          ['set-test']\n"
"  --model ['perm' | 'mperm='<value>] ['perm-count']\n"
"          [{fisher | fisher-midp | trend-only}] ['set-test']\n"
"          [{dom | rec | gen | trend}]\n"
"    Basic association analysis report.\n"
"    Given a case/control phenotype, --assoc performs a 1df chi-square allelic\n"
"    test, while --model performs 4 other tests as well (1df dominant gene\n"
"    action, 1df recessive gene action, 2df genotypic, Cochran-Armitage trend).\n"
"    * With 'fisher'/'fisher-midp', Fisher's exact test is used to generate\n"
"      p-values.  'fisher-midp' also applies Lancaster's mid-p adjustment.\n"
"    * 'perm' causes an adaptive permutation test to be performed.\n"
"    * 'mperm='<value> causes a max(T) permutation test with the specified\n"
"      number of replications to be performed.\n"
	       /*
"    * 'genedrop' causes offspring genotypes to be regenerated via gene-dropping\n"
"      in the permutation test.\n"
	       */
"    * 'perm-count' causes the permutation test report to include counts instead\n"
"      of frequencies.\n"
"    * 'counts' causes --assoc to report allele counts instead of frequencies.\n"
"    * 'set-test' tests the significance of variant sets.  Requires permutation;\n"
"      can be customized with --set-p/--set-r2/--set-max.\n"
"    * 'dom', 'rec', 'gen', and 'trend' force the corresponding test to be used\n"
"      as the basis for --model permutation.  (By default, the most significant\n"
"      result among the allelic, dominant, and recessive tests is used.)\n"
"    * 'trend-only' causes only the trend test to be performed.\n"
"    Given a quantitative phenotype, --assoc normally performs a Wald test.\n"
"    * In this case, the 'qt-means' modifier causes trait means and standard\n"
"      deviations stratified by genotype to be reported as well.\n"
"    * 'lin' causes the Lin statistic to be computed, and makes it the basis for\n"
"      multiple-testing corrections and permutation tests.\n"
"    Several other flags (most notably, --aperm) can be used to customize the\n"
"    permutation test.\n\n"
	       );
    help_print("mh\tbd\tmh2\thomog\tcmh", &help_ctrl, 1,
"  --mh ['perm' | 'mperm='<value>] ['perm-count'] ['set-test']\n"
"    (alias: --cmh)\n"
"  --bd ['perm' | 'perm-bd' | 'mperm='<value>] ['perm-count'] ['set-test']\n"
"  --mh2\n"
"  --homog\n"
"    Given a case/control phenotype and a set of clusters, --mh computes 2x2xK\n"
"    Cochran-Mantel-Haenszel statistics for each variant, while --bd also\n"
"    performs the Breslow-Day test for odds ratio homogeneity.  Permutation and\n"
"    variant set testing based on the CMH (default) or Breslow-Day (when\n"
"    'perm-bd' is present) statistic are supported.\n"
"    The following similar analyses are also available:\n"
"    * --mh2 swaps the roles of case/control status and cluster membership,\n"
"      performing a phenotype-stratified IxJxK Cochran-Mantel-Haenszel test on\n"
"      association between cluster assignments and genotypes.\n"
"    * --homog executes an alternative to the Breslow-Day test, based on\n"
"      partitioning of the chi-square statistic.\n\n"
	       );
    help_print("gxe\tmcovar", &help_ctrl, 1,
"  --gxe [covariate index]\n"
"    Given both a quantitative phenotype and a case/control covariate loaded\n"
"    with --covar defining two groups, --gxe compares the regression coefficient\n"
"    derived from considering only members of one group to the regression\n"
"    coefficient derived from considering only members of the other.  By\n"
"    default, the first covariate in the --covar file defines the groups; use\n"
"    e.g. \"--gxe 3\" to base them on the third covariate instead.\n\n"
	       );
    help_print("linear\tlogistic\tperm\tmperm\tperm-count\tset-test\tgenotypic\thethom\tdominant\trecessive\tno-snp\thide-covar\tsex\tno-x-sex\tinteraction\tstandard-beta\tbeta", &help_ctrl, 1,
#ifndef NOLAPACK
"  --linear ['perm' | 'mperm='<value>] ['perm-count'] ['set-test']\n"
"           [{genotypic | hethom | dominant | recessive | no-snp}]\n"
"           ['hide-covar'] [{sex | no-x-sex}] ['interaction'] ['beta']\n"
"           ['standard-beta'] ['intercept']\n"
#endif
"  --logistic ['perm' | 'mperm='<value>] ['perm-count'] ['set-test']\n"
"             [{genotypic | hethom | dominant | recessive | no-snp}]\n"
"             ['hide-covar'] [{sex | no-x-sex}] ['interaction'] ['beta']\n"
"             ['intercept']\n"
"    Multi-covariate association analysis on a quantitative (--linear) or\n"
"    case/control (--logistic) phenotype.  Normally used with --covar.\n"
"    * 'perm' normally causes an adaptive permutation test to be performed on\n"
"      the main effect, while 'mperm='<value> starts a max(T) permutation test.\n"
"    * 'perm-count' causes the permutation test report to include counts instead\n"
"      of frequencies.\n"
"    * 'set-test' tests the significance of variant sets.  Requires permutation;\n"
"      can be customized with --set-p/--set-r2/--set-max.\n"
"    * The 'genotypic' modifier adds an additive effect/dominance deviation 2df\n"
"      joint test (0/1/2 and 0/1/0 coding), while 'hethom' uses 0/0/1 and 0/1/0\n"
"      coding instead.  If permutation is also requested, these modifiers cause\n"
"      permutation to be based on the joint test.\n"
"    * 'dominant' and 'recessive' specify a model assuming full dominance or\n"
"      recessiveness, respectively, for the A1 allele.\n"
"    * 'no-snp' causes regression to be performed only on the phenotype and the\n"
"      covariates, without reference to genomic data.  If permutation is also\n"
"      requested, results are reported for all covariates.\n"
"    * 'hide-covar' removes covariate-specific lines from the report.\n"
"    * By default, sex (male = 1, female = 0) is automatically added as a\n"
"      covariate on X chromosome variants, and nowhere else.  The 'sex' modifier\n"
"      causes it to be added everywhere, while 'no-x-sex' excludes it.\n"
"    * 'interaction' adds genotype x covariate interactions to the model.  This\n"
"      cannot be used with the usual permutation tests; use --tests to define\n"
"      the permutation test statistic instead.\n"
"    * 'intercept' causes intercepts to be included in the main report.\n"
"    * For logistic regressions, the 'beta' modifier causes regression\n"
"      coefficients instead of odds ratios to be reported.\n"
"    * With --linear, the 'standard-beta' modifier standardizes the phenotype\n"
"      and all predictors to zero mean and unit variance before regression.\n\n"
	       );
    help_print("dosage\twrite-dosage", &help_ctrl, 1,
"  --dosage <allele dosage file> ['noheader'] ['skip0='<i>] ['skip1='<j>]\n"
"           ['skip2='<k>] ['dose1'] ['format='<m>] ['Zout']\n"
"           [{occur | standard-beta}] ['sex'] ['case-control-freqs']\n"
"  --dosage <list file> list [{sepheader | noheader}] ['skip0='<i>]\n"
"           ['skip1='<j>] ['skip2='<k>] ['dose1'] ['format='<m>] ['Zout']\n"
"           [{occur | standard-beta}] ['sex'] ['case-control-freqs']\n"
"  --write-dosage\n"
"    Process (possibly gzipped) text files with variant-major allelic dosage\n"
"    data.  This cannot be used with a regular input fileset; instead, you must\n"
"    *only* specify a .fam and possibly a .map file, and you can't specify any\n"
"    other commands.\n"
"    * PLINK 2.0 will have first-class support for genotype probabilities.  An\n"
"      equivalent data import flag will be provided then, and --dosage will be\n"
"      retired.\n"
"    * By default, --dosage assumes that only one allelic dosage file should be\n"
"      loaded.  To specify multiple files,\n"
"      1. create a master list with one entry per line.  There are normally two\n"
"         supported formats for this list: just a filename per line, or variant\n"
"         batch numbers in the first column and filenames in the second.\n"
"      2. Provide the name of that list as the first --dosage parameter.\n"
"      3. Add the 'list' modifier.\n"
"    * By default, --dosage assumes the allelic dosage file(s) contain a header\n"
"      line, which has 'SNP' in column i+1, 'A1' in column i+j+2, 'A2' in column\n"
"      i+j+3, and sample FID/IIDs starting from column i+j+k+4.  (i/j/k are\n"
"      normally zero, but can be changed with 'skip0', 'skip1', and 'skip2'\n"
"      respectively.)  If such a header line is not present,\n"
"      * when all samples appear in the same order as they do in the .fam file,\n"
"        you can use the 'noheader' modifier.\n"
"      * Otherwise, use the 'sepheader' modifier, and append sample ID filenames\n"
"        to your 'list' file entries.\n"
"    * The 'format=' modifier lets you specify the number of values used to\n"
"      represent each dosage.  'format=1' normally indicates a single 0..2 A1\n"
"      expected count; 'dose1' modifies this to a 0..1 frequency.  'format=2'\n"
"      (the default) indicates a 0..1 homozygous A1 likelihood followed by a\n"
"      0..1 het likelihood, while 'format=3' indicates 0..1 hom A1, 0..1 het,\n"
"      0..1 hom A2.\n"
"    * 'Zout' causes the output file to be gzipped.\n"
"    * Normally, an association analysis is performed.  'standard-beta' and\n"
"      'sex' behave as they are supposed to with --linear/--logistic.\n"
"      'case-control-freqs' causes case and control allele frequencies to be\n"
"      reported separately.\n"
"    * There are three alternate modes which cause the association analysis to\n"
"      be skipped.\n"
"      * 'occur' requests a simple variant occurrence report.\n"
"      * --write-dosage causes a simple merged file matching the 'format'\n"
"        specification (not including 'dose1') to be generated.\n"
"      * --score applies a linear scoring system to the dosages.\n\n"
	       );
    help_print("lasso", &help_ctrl, 1,
"  --lasso <h2 estimate> [min lambda] ['report-zeroes']\n"
"    Estimate variant effect sizes via LASSO regression.  You must provide an\n"
"    additive heritability estimate to calibrate the regression.\n"
"    Note that this method may require a very large sample size (e.g. hundreds\n"
"    of thousands) to be effective on complex polygenic traits.\n\n"
	       );
    help_print("test-missing\tmissing\tperm\tmperm", &help_ctrl, 1,
"  --test-missing ['perm' | 'mperm='<value>] ['perm-count'] ['midp']\n"
"    Check for association between missingness and case/control status, using\n"
"    Fisher's exact test.  (Het. haploids are treated as missing.)\n"
"    The 'midp' modifier causes Lancaster's mid-p adjustment to be applied.\n\n"
	       );
    help_print("make-perm-pheno", &help_ctrl, 1,
"  --make-perm-pheno <ct>\n"
"    Generate phenotype permutations and write them to disk, without invoking an\n"
"    association test.\n\n"
	       );
#ifndef STABLE_BUILD
#ifndef NOLAPACK
    help_print("unrelated-heritability", &help_ctrl, 1,
"  --unrelated-heritability ['strict'] [tol] [initial covg] [initial covr]\n"
"    REML estimate of additive heritability, iterating with an accelerated\n"
"    variant of the EM algorithm until the rate of change of the log likelihood\n"
"    function is less than tol.  A scalar phenotype is required.\n"
"    * The 'strict' modifier forces regular EM to be used.  tol defaults to\n"
"      10^{-7}, genomic covariance prior defaults to 0.45, and residual\n"
"      covariance prior defaults to (1 - covg).\n"
"    * You can combine this with --grm/--grm-bin to reuse a previously computed\n"
"      relationship matrix.  (Specify the phenotype with --pheno.)\n"
"    * For more details about the method, see Vattikuti S, Guo J, Chow CC (2012)\n"
"      Heritability and Genetic Correlations Explained by Common SNPs for\n"
"      Metabolic Syndrome Traits.\n\n"
	       );
#endif
#endif
    help_print("tdt\tpoo\tperm\tmperm\tparentdt1\tparentdt2\tpat\tmat\tset-test", &help_ctrl, 1,
"  --tdt [{exact | exact-midp | poo}] ['perm' | 'mperm='<value>] ['perm-count']\n"
"        [{parentdt1 | parentdt2 | pat | mat}] ['set-test']\n"
"    Report transmission disequilibrium test statistics, given case/control\n"
"    phenotypes and pedigree information.\n"
"    * A Mendel error check is performed before the main tests; offending\n"
"      genotypes are treated as missing by this analysis.\n"
"    * By default, the basic TDT p-value is based on a chi-square test unless\n"
"      you request the exact binomial test with 'exact' or 'exact-midp'.\n"
"    * 'perm'/'mperm=' requests a family-based adaptive or max(T) permutation\n"
"      test.  By default, the permutation test statistic is the basic TDT\n"
"      p-value; 'parentdt1'/'parentdt2' cause parenTDT or combined test\n"
"      p-values, respectively, to be considered instead.\n"
"    * 'set-test' tests the significance of variant sets.  This cannot be used\n"
"      with exact tests for now.\n"
"    The 'poo' modifier causes a parent-of-origin analysis to be performed\n"
"    instead, with transmissions from heterozygous fathers and heterozygous\n"
"    mothers considered separately.\n"
"    * The parent-of-origin analysis does not currently support exact tests.\n"
"    * By default, the permutation test statistic is the absolute\n"
"      parent-of-origin test Z score; 'pat'/'mat' cause paternal or maternal TDT\n"
"      chi-square statistics, respectively, to be considered instead.\n\n"
	       );
#ifndef STABLE_BUILD
    help_print("dfam", &help_ctrl, 1,
"  --dfam ['no-unrelateds'] ['perm' | 'mperm='<value>] ['perm-count']\n"
"         ['set-test']\n"
"    Sib-TDT-based association test.  By default, clusters of unrelated\n"
"    individuals are included in the test; the 'no-unrelateds' modifier removes\n"
"    this component, leaving the original sib-TDT.\n\n"
	       );
#endif
    help_print("qfam\tqfam-between\tqfam-parents\tqfam-total", &help_ctrl, 1,
"  --qfam ['perm' | 'mperm='<value>] ['perm-count'] ['emp-se']\n"
"  --qfam-parents ['perm' | 'mperm='<value>] ['perm-count'] ['emp-se']\n"
"  --qfam-between ['perm' | 'mperm='<value>] ['perm-count'] ['emp-se']\n"
"  --qfam-total ['perm' | 'mperm='<value>] ['perm-count'] ['emp-se']\n"
"    QFAM family-based association test for quantitative traits.\n"
"    * A Mendel error check is performed before the main tests; offending\n"
"      genotypes are treated as missing by this analysis.\n"
"    * This procedure requires permutation.  'perm' and 'perm-count' have the\n"
"      usual meanings.  However, 'mperm='<value> just specifies a fixed number\n"
"      of permutations; the method does not support a proper max(T) test.\n"
"    * The 'emp-se' modifier adds BETA and EMP_SE (empirical standard error for\n"
"      beta) fields to the .perm output file.\n\n"
	       );
#ifndef STABLE_BUILD
    help_print("tucc", &help_ctrl, 1,
"  --tucc ['write-bed']\n"
"    By default, this generates a .tucc.ped file where, for each child with two\n"
"    parents in the dataset, one case sample is created with all the child's\n"
"    alleles, and one pseudocontrol sample is created with all untransmitted\n"
"    alleles.  (Both genotypes are set to missing when there's a Mendel error.)\n"
"    Add the 'write-bed' modifier to create a .tucc.bed + .tucc.bim + .tucc.fam\n"
"    fileset instead.\n\n"
	       );
#endif
    help_print("annotate", &help_ctrl, 1,
"  --annotate <PLINK report> ['attrib='<file>] ['ranges='<file>]\n"
"             ['filter='<file>] ['snps='<file>] [{NA | prune}] ['block']\n"
"             ['subset='<file>] ['minimal'] ['distance']\n"
"    Add annotations to a variant-based PLINK report.  This requires an\n"
"    annotation source:\n"
"    * 'attrib='<file> specifies a (possibly gzipped) attribute file.\n"
"    * 'ranges='<file> specifies a gene/range list file.\n"
"    (Both source types can be specified simultaneously.)  The following options\n"
"    are also supported:\n"
"    * 'filter='<file> causes only variants within one of the ranges in the file\n"
"      to be included in the new report.\n"
"    * 'snps='<file> causes only variants named in the file to be included in\n"
"      the new report.\n"
"    * The 'NA' modifier causes unannotated variants to have 'NA' instead of '.'\n"
"      in the new report's ANNOT column, while the 'prune' modifier excludes\n"
"      them entirely.\n"
"    * The 'block' modifier replaces the single ANNOT column with a 0/1-coded\n"
"      column for each possible annotation.\n"
"    * With 'ranges',\n"
"      * 'subset='<file> causes only intervals named in the subset file to be\n"
"        loaded from the ranges file.\n"
"      * interval annotations normally come with a parenthesized signed distance\n"
"        to the interval boundary (0 if the variant is located inside the\n"
"        interval; this is always true without --border).  They can be excluded\n"
"        with the 'minimal' modifier.\n"
"      * the 'distance' modifier adds 'DIST' and 'SGN' columns describing signed\n"
"        distance to the nearest interval.\n"
"    * When --pfilter is present, high p-values are filtered out.\n\n"
	       );
    help_print("clump", &help_ctrl, 1,
"  --clump <PLINK report filename(s)...>\n"
"    Process association analysis report(s) with 'SNP' and p-value columns,\n"
"    organizing results by LD-based clumps.  Multiple filenames can be separated\n"
"    by spaces or commas.\n\n"
	       );
    help_print("gene-report\tgene-list", &help_ctrl, 1,
"  --gene-report <PLINK report> <gene range file>\n"
"    Generate a gene-based report from a variant-based report.\n"
"    * When --pfilter is present, high p-values are filtered out.\n"
"    * When --extract (without 'range') is present, only variants named in the\n"
"      --extract file are considered.\n\n"
	       );
    help_print("meta-analysis", &help_ctrl, 1,
"  --meta-analysis <PLINK report filenames...>\n"
"  --meta-analysis <PLINK report filenames...> + [{logscale | qt}]\n"
"                  [{no-map | no-allele}] ['study'] ['report-all']\n"
"                  ['weighted-z']\n"
"    Perform a meta-analysis on several variant-based reports with 'SNP' and\n"
"    'SE' fields.\n"
"    * Normally, an 'OR' odds ratio field must also be present in each input\n"
"      file.  With 'logscale', 'BETA' log-odds values/regression coefficients\n"
"      are expected instead, but the generated report will still contain odds\n"
"      ratio estimates.  With 'qt', both input and output values are regression\n"
"      betas.\n"
"    * 'CHR', 'BP', and 'A1' fields are also normally required.  'no-map' causes\n"
"      them to all be ignored, while 'no-allele' causes just 'A1' to be ignored.\n"
"    * If 'A2' fields are present, and neither 'no-map' nor 'no-allele' was\n"
"      specified, A1/A2 allele flips are handled properly.  Otherwise, A1\n"
"      mismatches are thrown out.\n"
"    * 'study' causes study-specific effect estimates to be collated in the\n"
"      meta-analysis report.\n"
"    * 'report-all' causes variants present in only a single input file to be\n"
"      included in the meta-analysis report.\n"
"    * 'weighted-z' requests weighted Z-score-based p-values (as computed by the\n"
"      Abecasis Lab's METAL software) in addition to the usual inverse\n"
"      variance-based analysis.  This requires P and effective sample size\n"
"      fields.\n"
"    * When --extract (without 'range') is present, only variants named in the\n"
"      --extract file are considered.\n"
"    * Unless 'no-map' is specified, chromosome filters are also respected.\n\n"
	       );
    help_print("fast-epistasis\tepistasis\tset-test\tset-by-all\tcase-only\tnop\tepistasis-summary-merge", &help_ctrl, 1,
"  --fast-epistasis [{boost | joint-effects | no-ueki}] ['case-only']\n"
"                   [{set-by-set | set-by-all}] ['nop']\n"
"  --epistasis [{set-by-set | set-by-all}]\n"
"    Scan for epistatic interactions.  --fast-epistasis inspects 3x3 joint\n"
"    genotype count tables and only applies to case/control phenotypes, while\n"
"    --epistasis performs linear or logistic regression.\n"
"    * By default, --fast-epistasis uses the PLINK 1.07 allele-based test.  Two\n"
"      newer tests are now supported: 'boost' invokes the likelihood ratio test\n"
"      introduced by Wan X et al. (2010) BOOST: A Fast Approach to Detecting\n"
"      Gene-Gene Interactions in Genome-wide Case-Control Studies, while\n"
"      'joint-effects' applies the joint effects test introduced in Ueki M,\n"
"      Cordell HJ (2012) Improved statistics for genome-wide interaction\n"
"      analysis.\n"
"    * The original --fast-epistasis test normally applies the variance and\n"
"      empty cell corrections suggested by Ueki and Cordell's paper.  To disable\n"
"      them, use the 'no-ueki' modifier.\n"
"    * 'case-only' requests a case-only instead of a case/control test.\n"
"    * By default, all pairs of variants across the entire genome are tested.\n"
"      To just test pairs of variants within a single set, add the 'set-by-set'\n"
"      modifier and load exactly one set with --set/--make-set; with exactly two\n"
"      sets loaded, all variants in one set are tested against all variants in\n"
"      the other.  'set-by-all' tests all variants in one set against the entire\n"
"      genome instead.\n"
"    * 'nop' strips p-values from the main report.\n"
"    * These computations can be subdivided with --parallel; however...\n"
"  --epistasis-summary-merge <common file prefix> <ct>\n"
"    When a --[fast-]epistasis job is subdivided with --parallel, the main\n"
"    report can be assembled at the end by applying Unix 'cat' in the usual\n"
"    manner, but the .summary.1, .summary.2, ... files may require a specialized\n"
"    merge.  --epistasis-summary-merge takes care of the latter.\n\n"
	       );
    help_print("twolocus", &help_ctrl, 1,
"  --twolocus <variant ID> <variant ID>\n"
"    Two-locus joint genotype count report.\n\n"
	       );
    help_print("score\tscore-no-mean-imputation", &help_ctrl, 1,
"  --score <filename> [i] [j] [k] ['header'] [{sum | no-sum}]\n"
"          [{no-mean-imputation | center}] ['include-cnt'] ['double-dosage']\n"
"    Apply a linear scoring system to each sample.\n"
"    The input file should have one line per scored variant.  Variant IDs are\n"
"    read from column #i, allele codes are read from column #j, and scores are\n"
"    read from column #k, where i defaults to 1, j defaults to i+1, and k\n"
"    defaults to j+1.\n"
"    * The 'header' modifier causes the first nonempty line of the input file to\n"
"      be ignored; otherwise, --score assumes there is no header line.\n"
"    * By default, final scores are averages of the valid per-variant scores.\n"
"      The 'sum' modifier causes sums to be reported instead.  (This cannot be\n"
"      used with 'no-mean-imputation'.  And for backward compatibility, 'sum' is\n"
"      automatically on with dosage data unless 'no-sum' is specified.)\n"
"    * By default, copies of the unnamed allele contribute zero to score, while\n"
"      missing genotypes contribute an amount proportional to the loaded (via\n"
"      --read-freq) or imputed allele frequency.  To throw out missing\n"
"      observations instead (decreasing the denominator in the final average\n"
"      when this happens), use the 'no-mean-imputation' modifier.\n"
"    * Alternatively, you can use the 'center' modifier to shift all scores to\n"
"      mean zero.\n"
"    * This command can be used with dosage data.  By default, the 'CNT' column\n"
"      is omitted from the output file in this case; use 'include-cnt' to keep\n"
"      it.  Also, note that scores are multiplied by 0..1 dosages, not 0..2\n"
"      diploid allele counts, unless the 'double-dosage' modifier is present.\n\n"
	       );
#if defined __cplusplus && !defined _WIN32
    help_print("R\tR-debug", &help_ctrl, 1,
"  --R <R script file> ['debug']\n"
"    Connect to a Rserve (preferably version 1.7 or later) background process,\n"
"    and execute the Rplink function defined in the input file.  (Unless the\n"
"    'debug' modifier is present; in that case, the R commands that PLINK would\n"
"    have tried to execute are logged to a file.)\n\n"
	       );
#endif
#ifndef STABLE_BUILD
    help_print("cnv-make-map", &help_ctrl, 1,
"  --cnv-make-map ['short']\n"
"    Given a .cnv file, this generates the corresponding .cnv.map file needed\n"
"    by " PROG_NAME_CAPS "'s other CNV analysis commands.  The 'short' modifier causes\n"
"    causes entries needed by old PLINK versions to be omitted.  (This is now\n"
"    automatically invoked, with 'short', when necessary.)\n\n"
	       );
    help_print("cnv-write", &help_ctrl, 1,
"  --cnv-write ['freq']\n"
"    Write a new .cnv fileset, after applying all requested filters.  The 'freq'\n"
"    modifier (which must be used with --cnv-freq-method2) causes an additional\n"
"    'FREQ' field to be written with CNV-CNV overlap counts.\n\n"
	       );
    /*
    help_print("cnv-check-no-overlap", &help_ctrl, 1,
"  --cnv-check-no-overlap\n"
"    Given a .cnv fileset, this checks for within-sample CNV overlaps.\n\n"
	       );
    */
    help_print("cnv-indiv-perm\tcnv-test\tcnv-test-region\tcnv-enrichment-test\tmperm\tcnv-test-1sided\tcnv-test-2sided", &help_ctrl, 1,
"  --cnv-indiv-perm <permutation count>\n"
"  --cnv-test [{1sided | 2sided}] <permutation count>\n"
"  --cnv-test-region <permutation count>\n"
"  --cnv-enrichment-test [permutation count]\n"
"    Given a .cnv fileset,\n"
"    * --cnv-indiv-perm performs a case/control CNV burden test.\n"
"    * --cnv-test performs a basic permutation-based CNV association test.  By\n"
"      default, this is 1-sided for case/control phenotypes and 2-sided for\n"
"      quantitative traits; you can use '1sided'/'2sided' to force the opposite.\n"
"    * --cnv-test-region reports case/ctrl association test results by region.\n"
"    * --cnv-enrichment-test performs Raychaudhuri et al.'s geneset enrichment\n"
"      test.  Gene locations must be loaded with --cnv-count.\n\n"
	       );
#endif
    help_print("write-var-ranges", &help_ctrl, 1,
"  --write-var-ranges <block ct>\n"
"    Divide the set of variants into equal-size blocks.  (Can be used with\n"
"    --snps to split a job across multiple machines.)\n\n"
	       );
    if (!param_ct) {
      fputs(
"The following other flags are supported.  (Order of operations is described at\n"
"https://www.cog-genomics.org/plink/1.9/order .)\n"
, stdout);
    }
    help_print("script\trerun", &help_ctrl, 0,
"  --script <fname> : Include command-line options from file.\n"
"  --rerun [log]    : Rerun commands in log (default '" PROG_NAME_STR ".log').\n"
	       );
    help_print("version", &help_ctrl, 0,
"  --version        : Display only version number before exiting.\n"
	       );
    help_print("silent\tgplink", &help_ctrl, 0,
"  --silent         : Suppress output to console.\n"
"  --gplink         : Reserved for interoperation with gPLINK.\n"
	       );
    help_print("missing-genotype", &help_ctrl, 0,
"  --missing-genotype <char> : Set missing genotype code (normally '0').\n"
	       );
    help_print("vcf\tbcf\tdouble-id\tconst-fid\tid-delim", &help_ctrl, 0,
"  --double-id          : Set both FIDs and IIDs to the VCF/BCF sample ID.\n"
"  --const-fid [ID]     : Set all FIDs to the given constant (default '0').\n"
"  --id-delim [d]       : Parse sample IDs as <FID><d><IID> (default delim '_').\n"
	       );
    help_print("vcf\tbcf\tid-delim\tvcf-idspace-to", &help_ctrl, 0,
"  --vcf-idspace-to <c> : Convert spaces in sample IDs to the given character.\n"
	       );
    help_print("vcf\tbcf\tbiallelic-only\tvcf-min-qual\tvcf-filter\tvcf-half-call\tvcf-min-gq\tvcf-min-gp\tvcf-require-gt", &help_ctrl, 0,
"  --biallelic-only ['strict'] ['list'] : Skip VCF variants with 2+ ALT alleles.\n"
"  --vcf-min-qual <val>           : Skip VCF variants with low/missing QUAL.\n"
"  --vcf-filter [exception(s)...] : Skip variants which have FILTER failures.\n"
"  --vcf-require-gt               : Skip variants with no GT field.\n"
"  --vcf-min-gq <val>             : No-call a genotype when GQ is below the\n"
"                                   given threshold.\n"
"  --vcf-min-gp <val>             : No-call a genotype when 0-1 scaled GP is\n"
"                                   below the given threshold.\n"
"  --vcf-half-call <m>  : Specify how '0/.' and similar VCF GT values should be\n"
"                         handled.  The following four modes are supported:\n"
"                         * 'error'/'e' (default) errors out and reports line #.\n"
"                         * 'haploid'/'h' treats them as haploid calls.\n"
"                         * 'missing'/'m' treats them as missing.\n"
"                         * 'reference'/'r' treats the missing value as 0.\n"
	       );
    help_print("oxford-single-chr\tdata\tgen", &help_ctrl, 0,
"  --oxford-single-chr <chr nm> : Specify single-chromosome .gen file with\n"
"                                 ignorable first column.\n"
	       );
    help_print("oxford-pheno-name\tdata\tsample", &help_ctrl, 0,
"  --oxford-pheno-name <col nm> : Import named phenotype from the .sample file.\n"
	       );
    help_print("hard-call-threshold\tmissing-code\tmissing_code\tdata\tgen\tbgen\tsample", &help_ctrl, 0,
"  --hard-call-threshold <val>  : When an Oxford-format fileset is loaded, calls\n"
"  --hard-call-threshold random   with uncertainty level greater than 0.1 are\n"
"                                 normally treated as missing.  You can adjust\n"
"                                 this threshold by providing a numeric\n"
"                                 parameter, or randomize all calls with\n"
"                                 'random'.\n"
"  --missing-code [string list] : Comma-delimited list of missing phenotype\n"
"    (alias: --missing_code)      values for Oxford-format filesets (def. 'NA').\n"
	       );
    help_print("simulate\tsimulate-ncases\tsimulate-ncontrols\tsimulate-prevalence", &help_ctrl, 0,
"  --simulate-ncases <num>   : Set --simulate case count (default 1000).\n"
"  --simulate-ncontrols <n>  : Set --simulate control count (default 1000).\n"
"  --simulate-prevalence <p> : Set --simulate disease prevalence (default 0.01).\n"
	       );
    help_print("simulate-qt\tsimulate-n", &help_ctrl, 0,
"  --simulate-n <num>        : Set --simulate-qt sample count (default 1000).\n"
	       );
    help_print("simulate\tsimulate-qt\tsimulate-label\tsimulate-missing", &help_ctrl, 0,
"  --simulate-label <prefix> : Set --simulate[-qt] FID/IID name prefix.\n"
"  --simulate-missing <freq> : Set --simulate[-qt] missing genotype frequency.\n"
	       );
    help_print("allow-extra-chr\taec", &help_ctrl, 0,
"  --allow-extra-chr ['0']   : Permit unrecognized chromosome codes.  The '0'\n"
"    (alias: --aec)            modifier causes them to be treated as if they had\n"
"                              been set to zero.\n"
               );
    help_print("chr-set\tcow\tdog\thorse\tmouse\trice\tsheep\tautosome-num", &help_ctrl, 0,
"  --chr-set <autosome ct> ['no-x'] ['no-y'] ['no-xy'] ['no-mt'] :\n"
"    Specify a nonhuman chromosome set.  The first parameter sets the number of\n"
"    diploid autosome pairs if positive, or haploid chromosomes if negative.\n"
"    Given diploid autosomes, the remaining modifiers indicate the absence of\n"
"    the named non-autosomal chromosomes.\n"
"  --cow/--dog/--horse/--mouse/--rice/--sheep : Shortcuts for those species.\n"
"  --autosome-num <value>    : Alias for \"--chr-set <value> no-y no-xy no-mt\".\n"
	       );
    help_print("cm-map\tzero-cms\tupdate-cm", &help_ctrl, 0,
"  --cm-map <fname pattern> [chr] : Use SHAPEIT-format recombination maps to set\n"
"                                   centimorgan positions.  To process more than\n"
"                                   one chromosome, include a '@' in the first\n"
"                                   parameter where the chrom. number belongs,\n"
"                                   e.g. 'genetic_map_chr@_combined_b37.txt'.\n"
"  --zero-cms         : Zero out centimorgan positions.\n"
	       );
    help_print("allow-no-samples\tallow-no-vars", &help_ctrl, 0,
"  --allow-no-samples : Allow the input fileset to contain no samples.\n"
"  --allow-no-vars    : Allow the input fileset to contain no variants.\n"
	       );
    help_print("pheno\tall-pheno\tmpheno\tpheno-name\tpheno-merge", &help_ctrl, 0,
"  --pheno <fname>  : Load phenotype data from the specified file, instead of\n"
"                     using the values in the main input fileset.\n"
"  --all-pheno      : For basic association tests, loop through all phenotypes\n"
"                     in --pheno file.\n"
"  --mpheno <n>     : Load phenotype from column (n+2) in --pheno file.\n"
"  --pheno-name <c> : If --pheno file has a header row, use column with the\n"
"                     given name.\n"
"  --pheno-merge    : When the main input fileset contains an phenotype value\n"
"                     for a sample, but the --pheno file does not, use the\n"
"                     original value instead of treating the phenotype as\n"
"                     missing.\n"
	       );
    help_print("missing-phenotype\t1", &help_ctrl, 0,
"  --missing-phenotype <v> : Set missing phenotype value (normally -9).\n"
"  --1                     : Expect case/control phenotypes to be coded as\n"
"                            0 = control, 1 = case, instead of the usual\n"
"                            0 = missing, 1 = control, 2 = case.  This also\n"
"                            forces phenotypes to be interpreted as case/ctrl.\n"
	       );
    help_print("make-pheno\tpheno", &help_ctrl, 0,
"  --make-pheno <fn> <val> : Define a new case/control phenotype.  If the val\n"
"                            parameter is '*', all samples listed in the given\n"
"                            file are cases, and everyone else is a control.\n"
"                            (Note that, in some shells, it is necessary to\n"
"                            surround the * with quotes.)\n"
"                            Otherwise, all samples with third column entry\n"
"                            equal to the val parameter are cases, and all other\n"
"                            samples mentioned in the file are controls.\n"
	       );
    help_print("tail-pheno\tgroupdist\tpheno", &help_ctrl, 0,
"  --tail-pheno <Lt> [Hbt] : Downcode a scalar phenotype to a case/control\n"
"                            phenotype.  All samples with phenotype values\n"
"                            greater than Hbt are cases, and all with values\n"
"                            less than or equal to Lt are controls.  If Hbt is\n"
"                            unspecified, it is equal to Lt; otherwise,\n"
"                            in-between phenotype values are set to missing.\n"
	       );
    help_print("covar\tcovar-name\tcovar-number\tno-const-covar\tallow-no-covars", &help_ctrl, 0,
"  --covar <filename> ['keep-pheno-on-missing-cov'] : Specify covariate file.\n"
"  --covar-name <...>       : Specify covariate(s) in --covar file by name.\n"
"                             Separate multiple names with spaces or commas, and\n"
"                             use dashes to designate ranges.\n"
"  --covar-number <...>     : Specify covariate(s) in --covar file by index.\n"
"  --no-const-covar         : Exclude constant covariates.\n"
"  --allow-no-covars        : Allow no covariates to be loaded from --covar\n"
"                             file.\n"
	       );
    help_print("within\tmwithin\tfamily", &help_ctrl, 0,
"  --within <f> ['keep-NA'] : Specify initial cluster assignments.\n"
"  --mwithin <n>            : Load cluster assignments from column n+2.\n"
"  --family                 : Create a cluster for each family ID.\n"
	       );
    help_print("loop-assoc", &help_ctrl, 0,
"  --loop-assoc <f> ['keep-NA']  : Run specified case/control association\n"
"                                  commands once for each cluster in the file,\n"
"                                  using cluster membership as the phenotype.\n"
	       );
    help_print("set\tset-names\tsubset\tset-collapse-all\tmake-set-collapse-all\tcomplement-sets\tmake-set-complement-all\tmake-set\tmake-set-border\tborder\tmake-set-collapse-group\t--make-set-complement-group", &help_ctrl, 0,
"  --set <filename>              : Load sets from a .set file.\n"
"  --set-names <name(s)...>      : Load only sets named on the command line.\n"
"                                  Use spaces to separate multiple names.\n"
"  --subset <filename>           : Load only sets named in the given text file.\n"
"  --set-collapse-all <set name> : Merge all sets.\n"
"  --complement-sets             : Invert all sets.  (Names gain 'C_' prefixes.)\n"
"  --make-set-complement-all <s> : --set-collapse-all + inversion.\n"
"  --make-set <filename>         : Define sets from a list of named bp ranges.\n"
"  --make-set-border <kbs>       : Stretch regions in --make-set file.\n"
"  --make-set-collapse-group     : Define sets from groups instead of sets in\n"
"                                  --make-set file.\n"
	       );
    help_print("keep\tremove\tkeep-fam\tremove-fam", &help_ctrl, 0,
"  --keep <filename>       : Exclude all samples not named in the file.\n"
"  --remove <filename>     : Exclude all samples named in the file.\n"
"  --keep-fam <filename>   : Exclude all families not named in the file.\n"
"  --remove-fam <filename> : Exclude all families named in the file.\n"
	       );
    help_print("extract\texclude\trange", &help_ctrl, 0,
"  --extract ['range'] <f> : Exclude all variants not named in the file.\n"
"  --exclude ['range'] <f> : Exclude all variants named in the file.\n"
	       );
    help_print("keep-clusters\tkeep-cluster-names\tremove-clusters\tremove-cluster-names", &help_ctrl, 0,
"  --keep-clusters <filename>          : These can be used individually or in\n"
"  --keep-cluster-names <name(s)...>     combination to define a list of\n"
"                                        clusters to keep; all samples not in a\n"
"                                        cluster in that list are then excluded.\n"
"                                        Use spaces to separate cluster names\n"
"                                        for --keep-cluster-names.\n"
"  --remove-clusters <filename>        : Exclude all clusters named in the file.\n"
"  --remove-cluster-names <name(s)...> : Exclude the named clusters.\n"
	       );
    help_print("set\tmake-set\tgene\tgene-all", &help_ctrl, 0,
"  --gene <sets...> : Exclude variants not in a set named on the command line.\n"
"                     (Separate multiple set names with spaces.)\n"
"  --gene-all       : Exclude variants which aren't a member of any set.  (PLINK\n"
"                     1.07 automatically did this under some circumstances.)\n"
	       );
    help_print("attrib\tattrib-indiv", &help_ctrl, 0,
"  --attrib <f> [att lst] : Given a file assigning attributes to variants, and a\n"
"  --attrib-indiv <f> [a]   comma-delimited list (with no whitespace) of\n"
"                           attribute names, remove variants/samples which are\n"
"                           either missing from the file or don't have any of\n"
"                           the listed attributes.  If some attribute names in\n"
"                           the list are preceded by '-', they are treated as\n"
"                           \"negative match conditions\" instead: variants with\n"
"                           at least one negative match attribute are removed.\n"
"                           The first character in the list cannot be a '-', due\n"
"                           to how command-line parsing works; add a comma in\n"
"                           front to get around this.\n"
	       );
    help_print("chr\tnot-chr\tchr-excl\tfrom-bp\tto-bp\tfrom-kb\tto-kb\tfrom-mb\tto-mb", &help_ctrl, 0,
"  --chr <chrs...>  : Exclude all variants not on the given chromosome(s).\n"
"                     Valid choices for humans are 0 (unplaced), 1-22, X, Y, XY,\n"
"                     and MT.  Separate multiple chromosomes with spaces and/or\n"
"                     commas, and use a dash (no adjacent spaces permitted) to\n"
"                     denote a range, e.g. \"--chr 1-4, 22, xy\".\n"
"  --not-chr <...>  : Reverse of --chr (exclude variants on listed chromosomes).\n"
	       );
    help_print("autosome\tautosome-xy\tchr\tnot-chr\tchr-excl", &help_ctrl, 0,
"  --autosome       : Exclude all non-autosomal variants.\n"
"  --autosome-xy    : Exclude all non-autosomal variants, except those with\n"
"                     chromosome code XY (pseudo-autosomal region of X).\n"
	       );
    help_print("snps-only", &help_ctrl, 0,
"  --snps-only ['just-acgt'] : Exclude non-SNP variants.  By default, SNP = both\n"
"                              allele codes are single-character; 'just-acgt'\n"
"                              restricts codes to {A,C,G,T,a,c,g,t,<missing>}.\n"
	       );
    help_print("from\tto\tsnp\twindow\tfrom-bp\tto-bp\tfrom-kb\tto-kb\tfrom-mb\tto-mb\texclude-snp\textract-snp", &help_ctrl, 0,
"  --from <var ID>  : Use ID(s) to specify a variant range to load.  When used\n"
"  --to   <var ID>    together, both variants must be on the same chromosome.\n"
"  --snp  <var ID>  : Specify a single variant to load.\n"
"  --exclude-snp <> : Specify a single variant to exclude.\n"
"  --window  <kbs>  : With --snp or --exclude-snp, loads/excludes all variants\n"
"                     within half the specified kb distance of the named one.\n"
"  --from-bp <pos>  : Use physical position(s) to define a variant range to\n"
"  --to-bp   <pos>    load.  --from-kb/--to-kb/--from-mb/--to-mb allow decimal\n"
"  --from-kb <pos>    values.  You must also specify a single chromosome (using\n"
"  --to-kb   <pos>    e.g. --chr) when using these flags.\n"
"  --from-mb <pos>\n"
"  --to-mb   <pos>\n"
	       );
    help_print("snps\texclude-snps", &help_ctrl, 0,
"  --snps <var IDs...>  : Use IDs to specify variant range(s) to load or\n"
"  --exclude-snps <...>   exclude.  E.g. \"--snps rs1111-rs2222, rs3333, rs4444\".\n"
	       );
    help_print("thin\tthin-count", &help_ctrl, 0,
"  --thin <p>       : Randomly remove variants, retaining each with prob. p.\n"
"  --thin-count <n> : Randomly remove variants until n of them remain.\n"
	       );
    help_print("bp-space\tthin", &help_ctrl, 0,
"  --bp-space <bps> : Remove variants so that each pair is no closer than the\n"
"                     given bp distance.  (Equivalent to VCFtools --thin.)\n"
	       );
    help_print("thin-indiv\tthin-indiv-count\tmax-indv", &help_ctrl, 0,
"  --thin-indiv <p>         : Randomly remove samples, retaining with prob. p.\n"
"  --thin-indiv-count <n>   : Randomly remove samples until n of them remain.\n"
	       );
    help_print("filter\tmfilter", &help_ctrl, 0,
"  --filter <f> <val(s)...> : Exclude all samples without a 3rd column entry in\n"
"                             the given file matching one of the given\n"
"                             space-separated value(s).\n"
"  --mfilter <n>            : Match against (n+2)th column instead.\n"
	       );
    help_print("geno\tmind\toblig-clusters\toblig-missing", &help_ctrl, 0,
"  --geno [val]     : Exclude variants with missing call frequencies greater\n"
"                     than a threshold (default 0.1).  (Note that the default\n"
"                     threshold is only applied if --geno is invoked without a\n"
"                     parameter; when --geno is not invoked, no per-variant\n"
"                     missing call frequency ceiling is enforced at all.  Other\n"
"                     inclusion/exclusion default thresholds work the same way.)\n"
"  --mind [val]     : Exclude samples with missing call frequencies greater than\n"
"                     a threshold (default 0.1).\n"
	       );
    help_print("oblig-clusters\toblig-missing", &help_ctrl, 0,
"  --oblig-missing <f1> <f2> : Specify blocks of missing genotype calls for\n"
"                              --geno/--mind to ignore.  The first file should\n"
"                              have variant IDs in the first column and block\n"
"                              IDs in the second, while the second file should\n"
"                              have FIDs in the first column, IIDs in the\n"
"                              second, and block IDs in the third.\n"
	       );
    help_print("prune", &help_ctrl, 0,
"  --prune             : Remove samples with missing phenotypes.\n"
	       );
    help_print("maf\tmax-maf\tmac\tmin-ac\tmax-mac\tmax-ac", &help_ctrl, 0,
"  --maf [freq]        : Exclude variants with minor allele frequency lower than\n"
"                        a threshold (default 0.01).\n"
"  --max-maf <freq>    : Exclude variants with MAF greater than the threshold.\n"
"  --mac <ct>          : Exclude variants with minor allele count lower than the\n"
"    (alias: --min-ac)   given threshold.\n"
"  --max-mac <ct>      : Exclude variants with minor allele count greater than\n"
"    (alias: --max-ac)   the given threshold.\n"
	       );
    help_print("maf-succ", &help_ctrl, 0,
"  --maf-succ       : Rule of succession MAF estimation (used in EIGENSOFT).\n"
"                     Given j observations of one allele and k >= j observations\n"
"                     of the other, infer a MAF of (j+1) / (j+k+2), rather than\n"
"                     the default j / (j+k).\n"
	       );
    help_print("read-freq\tupdate-freq", &help_ctrl, 0,
"  --read-freq <fn> : Estimate MAFs and heterozygote frequencies from the given\n"
"                     --freq[x] report, instead of the input fileset.\n"
	       );
    help_print("hwe\thwe-all\thwe2", &help_ctrl, 0,
"  --hwe <p> ['midp'] ['include-nonctrl'] : Exclude variants with Hardy-Weinberg\n"
"                                           equilibrium exact test p-values\n"
"                                           below a threshold.\n"
	       );
    help_print("me\tme-exclude-one", &help_ctrl, 0,
"  --me <t> <v> ['var-first'] : Filter out trios and variants with Mendel error\n"
"                               rates exceeding the given thresholds.\n"
"  --me-exclude-one [ratio]   : Make --me exclude only one sample per trio.\n"
	       );
    help_print("qual-scores\tqual-threshold\tqual-max-threshold", &help_ctrl, 0,
"  --qual-scores <f> [qcol] [IDcol] [skip] : Filter out variants with\n"
"                                            out-of-range quality scores.\n"
"                                            Default range is now [0, \\infty ).\n"
"  --qual-threshold <min qual score>       : Set --qual-scores range floor.\n"
"  --qual-max-threshold <max qual score>   : Set --qual-scores range ceiling.\n"
	       );
    help_print("allow-no-sex\tmust-have-sex", &help_ctrl, 0,
"  --allow-no-sex   : Do not treat ambiguous-sex samples as having missing\n"
"                     phenotypes in analysis commands.  (Automatic /w --no-sex.)\n"
"  --must-have-sex  : Force ambiguous-sex phenotypes to missing on\n"
"                     --make-bed/--make-just-fam/--recode/--write-covar.\n"
	       );
    help_print("filter-cases\tfilter-controls", &help_ctrl, 0,
"  --filter-cases       : Include only cases in the current analysis.\n"
"  --filter-controls    : Include only controls.\n"
	       );
    help_print("filter-males\tfilter-females", &help_ctrl, 0,
"  --filter-males       : Include only males.\n"
"  --filter-females     : Include only females.\n"
	       );
    help_print("filter-founders\tfilter-nonfounders", &help_ctrl, 0,
"  --filter-founders    : Include only founders.\n"
"  --filter-nonfounders : Include only nonfounders.\n"
	       );
    help_print("nonfounders", &help_ctrl, 0,
"  --nonfounders        : Include nonfounders in allele freq/HWE calculations.\n"
	       );
    help_print("make-founders", &help_ctrl, 0,
"  --make-founders ['require-2-missing'] ['first'] :\n"
"    Clear parental IDs for those with 1+ missing parent(s).\n"
	       );
    help_print("recode\trecode-allele\texport", &help_ctrl, 0,
"  --recode-allele <fn> : With --recode A/A-transpose/AD, count alleles named in\n"
"                         the file (otherwise A1 alleles are always counted).\n"
	       );
    help_print("output-chr\tchr-output", &help_ctrl, 0,
"  --output-chr <MT code> : Set chromosome coding scheme in output files by\n"
"                           providing the desired human mitochondrial code.\n"
"                           (Options are '26', 'M', 'MT', '0M', 'chr26', 'chrM',\n"
"                           and 'chrMT'.)\n"
	       );
    help_print("output-missing-genotype\toutput-missing-phenotype", &help_ctrl, 0,
"  --output-missing-genotype <ch> : Set the code used to represent missing\n"
"                                   genotypes in output files (normally the\n"
"                                   --missing-genotype value).\n"
"  --output-missing-phenotype <s> : Set the string used to represent missing\n"
"                                   phenotypes in output files (normally the\n"
"                                   --missing-phenotype value).\n"
	       );
    help_print("zero-cluster", &help_ctrl, 0,
"  --zero-cluster <f> : In combination with --within/--family, set blocks of\n"
"                       genotype calls to missing.  The input file should have\n"
"                       variant IDs in the first column and cluster IDs in the\n"
"                       second.  This must now be used with --make-bed and no\n"
"                       other output commands.\n"
	       );
    help_print("set-hh-missing\tset-mixed-mt-missing", &help_ctrl, 0,
"  --set-hh-missing       : Cause --make-bed and --recode to set heterozygous\n"
"                           haploid genotypes to missing.\n"
"  --set-mixed-mt-missing : Cause --make-bed and --recode to set mixed MT\n"
"                           genotypes to missing.\n"
	       );
    help_print("split-x\tmerge-x\tset-hh-missing\t23file-convert-xy\t23file-make-xylist\tcheck-sex\timpute-sex", &help_ctrl, 0,
"  --split-x <bp1> <bp2> ['no-fail']\n"
"  --split-x <build> ['no-fail'] :\n"
"    Changes chromosome code of all chrX variants with bp position <= bp1 or >=\n"
"    bp2 to XY.  The following build codes are supported as shorthand:\n"
"    * 'b36'/'hg18' = NCBI 36, 2709521/154584237\n"
"    * 'b37'/'hg19' = GRCh37, 2699520/154931044\n"
"    * 'b38'/'hg38' = GRCh38, 2781479/155701383\n"
"    By default, PLINK errors out when no variants would be affected by\n"
"    --split-x; the 'no-fail' modifier (useful in scripts) overrides this.\n"
"  --merge-x ['no-fail'] : Merge XY chromosome back with X.\n"
	       );
    help_print("set-me-missing", &help_ctrl, 0,
"  --set-me-missing      : Cause --make-bed to set Mendel errors to missing.\n"
	       );
    help_print("fill-missing-a2", &help_ctrl, 0,
"  --fill-missing-a2     : Cause --make-bed to replace all missing calls with\n"
"                          homozygous A2 calls.\n"
	       );
    help_print("set-missing-snp-ids\tset-missing-nonsnp-ids\tset-missing-var-ids\tnew-id-max-allele-len\tmissing-var-code", &help_ctrl, 0,
"  --set-missing-var-ids <t>   : Given a template string with a '@' where the\n"
"                                chromosome code should go and '#' where the bp\n"
"                                coordinate belongs, --set-missing-var-ids\n"
"                                assigns chromosome-and-bp-based IDs to unnamed\n"
"                                variants.\n"
"                                You may also use '$1' and '$2' to refer to\n"
"                                allele names in the template string, and in\n"
"                                fact this becomes essential when multiple\n"
"                                variants share the same coordinate.\n"
"  --new-id-max-allele-len <n> : Specify maximum number of leading characters\n"
"                                from allele names to include in new variant IDs\n"
"                                (default 23).\n"
"  --missing-var-code <string> : Change unnamed variant code (default '.').\n"
	       );
    help_print("update-chr\tupdate-cm\tupdate-map\tupdate-name", &help_ctrl, 0,
"  --update-chr  <f> [chrcol] [IDcol]  [skip] : Update variant chromosome codes.\n"
"  --update-cm   <f> [cmcol]  [IDcol]  [skip] : Update centimorgan positions.\n"
"  --update-map  <f> [bpcol]  [IDcol]  [skip] : Update variant bp positions.\n"
"  --update-name <f> [newcol] [oldcol] [skip] : Update variant IDs.\n"
	       );
    help_print("update-alleles", &help_ctrl, 0,
"  --update-alleles <fname> : Update variant allele codes.\n"
	       );
    help_print("allele1234\talleleACGT\talleleacgt", &help_ctrl, 0,
"  --allele1234 ['multichar'] : Interpret/recode A/C/G/T alleles as 1/2/3/4.\n"
"                               With 'multichar', converts all A/C/G/Ts in\n"
"                               allele names to 1/2/3/4s.\n"
"  --alleleACGT ['multichar'] : Reverse of --allele1234.\n"
	       );
    help_print("update-ids\tupdate-parents\tupdate-sex\timpute-sex", &help_ctrl, 0,
"  --update-ids <f>     : Update sample IDs.\n"
"  --update-parents <f> : Update parental IDs.\n"
"  --update-sex <f> [n] : Update sexes.  Sex (1 or M = male, 2 or F = female, 0\n"
"                         = missing) is loaded from column n+2 (default n is 1).\n"
	       );
    help_print("flip\tflip-subset", &help_ctrl, 0,
"  --flip <filename>    : Flip alleles (A<->T, C<->G) for SNP IDs in the file.\n"
"  --flip-subset <fn>   : Only apply --flip to samples in --flip-subset file.\n"
	       );
    help_print("flip-scan\tflip-scan-window\tflip-scan-window-kb\tflip-scan-threshold\tld-window\tld-window-kb\tflipscan\tflipscan-window\tflipscan-window-kb\tflipscan-threshold", &help_ctrl, 0,
"  --flip-scan-window <ct+1> : Set --flip-scan max variant ct dist. (def. 10).\n"
"  --flip-scan-window-kb <x> : Set --flip-scan max kb distance (default 1000).\n"
"  --flip-scan-threshold <x> : Set --flip-scan min correlation (default 0.5).\n"
	       );
    help_print("keep-allele-order\treal-ref-alleles\tmake-bed\tmerge\tbmerge\tmerge-list\trecode", &help_ctrl, 0,
"  --keep-allele-order  : Keep the allele order defined in the .bim file,\n"
"  --real-ref-alleles     instead of forcing A2 to be the major allele.\n"
"                         --real-ref-alleles also removes 'PR' from the INFO\n"
"                         values emitted by --recode vcf{,-fid,-iid}.\n"
	       );
    help_print("a1-allele\treference-allele\tupdate-ref-allele\ta2-allele", &help_ctrl, 0,
"  --a1-allele <f> [a1col] [IDcol] [skip] : Force alleles in the file to A1.\n"
"  --a2-allele <filename> [a2col] [IDcol] [skip] :\n"
"    Force alleles in the file to A2.  (\"--a2-allele <VCF filename> 4 3 '#'\",\n"
"    which scrapes reference allele assignments from a VCF file, is especially\n"
"    useful.)\n"
	       );
    help_print("indiv-sort\tmerge\tbmerge\tmerge-list", &help_ctrl, 0,
"  --indiv-sort <m> [f] : Specify FID/IID sort order.  The following four modes\n"
"                         are supported:\n"
"                         * 'none'/'0' keeps samples in the order they were\n"
"                           loaded.  Default for non-merge operations.\n"
"                         * 'natural'/'n' invokes 'natural sort', e.g.\n"
"                           'id2' < 'ID3' < 'id10'.  Default when merging.\n"
"                         * 'ascii'/'a' sorts in ASCII order, e.g.\n"
"                           'ID3' < 'id10' < 'id2'.\n"
"                         * 'file'/'f' uses the order in the given file (named\n"
"                           in the second parameter).\n"
"                         For now, only --merge/--bmerge/--merge-list and\n"
"                         --make-bed/--make-just-fam respect this flag.\n"
	       );
    help_print("with-phenotype\tdummy-coding\twrite-covar", &help_ctrl, 0,
"  --with-phenotype ['no-parents'] [{no-sex | female-2}] :\n"
"    Include more sample info in new .cov file.\n"
"  --dummy-coding [N] ['no-round'] : Split categorical variables (n categories,\n"
"                                    2 < n <= N, default N is 49) into n-1\n"
"                                    binary dummy variables when writing\n"
"                                    covariate file.\n"
	       );
    help_print("merge\tbmerge\tmerge-list\tmerge-mode", &help_ctrl, 0,
"  --merge-mode <n>   : Adjust --[b]merge/--merge-list behavior based on a\n"
"                       numeric code.\n"
"                       1 (default) = ignore missing calls, otherwise difference\n"
"                                     -> missing\n"
"                       2 = only overwrite originally missing calls\n"
"                       3 = only overwrite when nonmissing in new file\n"
"                       4/5 = never overwrite and always overwrite, respectively\n"
"                       6 = report all mismatching calls without merging\n"
"                       7 = report mismatching nonmissing calls without merging\n"
	       );
    help_print("merge\tbmerge\tmerge-list\tmerge-mode\tmerge-equal-pos", &help_ctrl, 0,
"  --merge-equal-pos  : With --merge/--bmerge/--merge-list, merge variants with\n"
"                       different names but identical positions.  (Exception:\n"
"                       same-position chromosome code 0 variants aren't merged.)\n"
	       );
    help_print("mendel-duos\tmendel-multigen\tme\tmendel\ttdt\tset-me-missing", &help_ctrl, 0,
"  --mendel-duos      : Make Mendel error checks consider samples with only one\n"
"                       parent in the dataset.\n"
	       );
    help_print("mendel-duos\tmendel-multigen\tme\tmendel\ttdt\tset-me-missing\ttdt", &help_ctrl, 0,
"  --mendel-multigen  : Make Mendel error checks consider (great-)grandparental\n"
"                       genotypes when parental genotype data is missing.\n"
	       );
    help_print("r\tr2\tld-window-kb\tld-window-cm\tld-window-r2\tld-window\tld-snp\tld-snps\tld-snp-list", &help_ctrl, 0,
"  --ld-window <ct+1> : Set --r/--r2 max variant ct pairwise distance (usu. 10).\n"
"  --ld-window-kb <x> : Set --r/--r2 max kb pairwise distance (usually 1000).\n"
"  --ld-window-cm <x> : Set --r/--r2 max centimorgan pairwise distance.\n"
"  --ld-window-r2 <x> : Set threshold for --r2 report inclusion (usually 0.2).\n"
"  --ld-snp <var ID>  : Set first variant in all --r/--r2 pairs.\n"
"  --ld-snps <vID...> : Restrict first --r/--r2 variant to the given ranges.\n"
"  --ld-snp-list <f>  : Restrict first --r/--r2 var. to those named in the file.\n"
	       );
    help_print("list-all\ttag-kb\ttag-r2\ttag-mode2\tshow-tags", &help_ctrl, 0,
"  --list-all         : Generate the 'all' mode report when using --show-tags in\n"
"                       file mode.\n"
"  --tag-kb <kbs>     : Set --show-tags max tag kb distance (default 250).\n"
"  --tag-r2 <val>     : Set --show-tags min tag r-squared (default 0.8)\n"
"  --tag-mode2        : Use two-column --show-tags (file mode) I/O format.\n"
	       );
    help_print("indep\tindep-pairwise\tr\tr2\tflip-scan\tflipscan\tshow-tags\tld-xchr", &help_ctrl, 0,
"  --ld-xchr <code>   : Set chrX model for --indep[-pairwise], --r/--r2,\n"
"                       --flip-scan, and --show-tags.\n"
"                       1 (default) = males coded 0/1, females 0/1/2 (A1 dosage)\n"
"                       2 = males coded 0/2\n"
"                       3 = males coded 0/2, but females given double weighting\n"
	       );
    help_print("blocks\tblocks-max-kb\tblocks-min-maf\tblocks-strong-lowci\tblocks-strong-highci\tblocks-recomb-highci\tblocks-inform-frac\tld-window-kb", &help_ctrl, 0,
"  --blocks-max-kb <kbs>      : Set --blocks maximum haploblock span (def. 200).\n"
"  --blocks-min-maf <cutoff>  : Adjust --blocks MAF minimum (default 0.05).\n"
"  --blocks-strong-lowci <x>  : Set --blocks \"strong LD\" CI thresholds (defaults\n"
"  --blocks-strong-highci <x>   0.70 and 0.98).\n"
"  --blocks-recomb-highci <x> : Set 'recombination' CI threshold (default 0.90).\n"
"  --blocks-inform-frac <x>   : Force haploblock <strong LD pairs>:<total\n"
"                               informative pairs> ratios to be larger than this\n"
"                               value (default 0.95).\n"
	       );
    help_print("distance-wts\tdistance-exp\texponent\tdistance", &help_ctrl, 0,
"  --distance-wts exp=<x>        : When computing genomic distances, assign each\n"
"                                  variant a weight of (2q(1-q))^{-x}, where q\n"
"                                  is the loaded or inferred MAF.\n"
	       );
#ifndef STABLE_BUILD
    help_print("distance-wts\tdistance\tmake-grm-gz\tmake-grm-bin", &help_ctrl, 0,
"  --distance-wts <f> ['noheader'] : When computing genomic distances, assign\n"
"                                    each variant the weight specified in the\n"
"                                    file.\n"
	       );
#endif
    help_print("read-dists\tload-dists\tibs-test\tgroupdist\tregress-distance\tcluster\tneighbour\tneighbor", &help_ctrl, 0,
"  --read-dists <dist file> [id file] : Load a triangular binary distance matrix\n"
"                                       instead of recalculating from scratch.\n"
	       );
    help_print("ppc-gap\tmin\tmax\tgenome\tZ-genome", &help_ctrl, 0,
"  --ppc-gap <val>    : Minimum number of base pairs, in thousands, between\n"
"                       informative pairs of markers used in --genome PPC test.\n"
"                       500 if unspecified.\n"
"  --min <cutoff>     : Specify minimum PI_HAT for inclusion in --genome report.\n"
"  --max <cutoff>     : Specify maximum PI_HAT for inclusion in --genome report.\n"
	       );
    help_print("homozyg\thomozyg-match\tpool-size", &help_ctrl, 0,
"  --homozyg-match <> : Set minimum concordance across jointly homozygous\n"
"                       variants for a pairwise allelic match to be declared.\n"
"  --pool-size <ct>   : Set minimum size of pools in \"--homozyg group\" report.\n"
	       );
    help_print("read-genome\tcluster\tneighbour\tneighbor", &help_ctrl, 0,
"  --read-genome <fn> : Load --genome report for --cluster/--neighbour, instead\n"
"                       of recalculating IBS and PPC test p-values from scratch.\n"
	       );
    help_print("ppc\tmc\tmcc\tK\tk\tibm\tcluster", &help_ctrl, 0,
"  --ppc <p-val>    : Specify minimum PPC test p-value within a cluster.\n"
"  --mc <max size>  : Specify maximum cluster size.\n"
"  --mcc <c1> <c2>  : Specify maximum case and control counts per cluster.\n"
"  --K <min count>  : Specify minimum cluster count.\n"
"  --ibm <val>      : Specify minimum identity-by-missingness.\n"
	       );
    help_print("match\tmatch-type\tqmatch\tqt\tcluster", &help_ctrl, 0,
"  --match <f> [mv] : Use covariate values to restrict clustering.  Without\n"
"                     --match-type, two samples can only be in the same cluster\n"
"                     if all covariates match.  The optional second parameter\n"
"                     specifies a covariate value to treat as missing.\n"
"  --match-type <f> : Refine interpretation of --match file.  The --match-type\n"
"                     file is expected to be a single line with as many entries\n"
"                     as the --match file has covariates; '0' entries specify\n"
"                     \"negative matches\" (i.e. samples with equal covariate\n"
"                     values cannot be in the same cluster), '1' entries specify\n"
"                     \"positive matches\" (default), and '-1' causes the\n"
"                     corresponding covariate to be ignored.\n"
"  --qmatch <f> [m] : Force all members of a cluster to have similar\n"
"  --qt <fname>       quantitative covariate values.  The --qmatch file contains\n"
"                     the covariate values, while the --qt file is a list of\n"
"                     nonnegative tolerances (and '-1's marking covariates to\n"
"                     skip).\n"
	       );
#ifndef NOLAPACK
    help_print("pca\tpca-cluster-names\tpca-clusters", &help_ctrl, 0,
"  --pca-cluster-names <...> : These can be used individually or in combination\n"
"  --pca-clusters <fname>      to define a list of clusters to use in the basic\n"
"                              --pca computation.  (--pca-cluster-names expects\n"
"                              a space-delimited sequence of cluster names,\n"
"                              while --pca-clusters expects a file with one\n"
"                              cluster name per line.)  All samples outside\n"
"                              those clusters will then be projected on to the\n"
"                              calculated PCs.\n"
	       );
    help_print("cluster\tmds-plot\tmds-cluster", &help_ctrl, 0,
"  --mds-plot <dims> ['by-cluster'] ['eigendecomp'] ['eigvals'] :\n"
"    Multidimensional scaling analysis.  Requires --cluster.\n"
	       );
#endif
    help_print("cell\tmodel", &help_ctrl, 0,
"  --cell <thresh>  : Skip some --model tests when a contingency table entry is\n"
"                     smaller than the given threshold.\n"
	       );
    help_print("linear\tlogistic\tcondition\tcondition-list\tparameters\ttests\ttest-all", &help_ctrl, 0,
"  --condition <var ID> [{dominant | recessive}] : Add one variant as a --linear\n"
"                                                  or --logistic covariate.\n"
"  --condition-list <f> [{dominant | recessive}] : Add variants named in the\n"
"                                                  file as --linear/--logistic\n"
"                                                  covariates.\n"
"  --parameters <...>  : Include only the given covariates/interactions in the\n"
"                        --linear/--logistic models, identified by a list of\n"
"                        1-based indices and/or ranges of them.\n"
"  --tests <all> [...] : Perform a (joint) test on the specified term(s) in the\n"
"                        --linear/--logistic model, identified by 1-based\n"
"                        indices and/or ranges of them.  If permutation was\n"
"                        requested, it is based on this test.\n"
"                        * Note that, when --parameters is also present, the\n"
"                          indices refer to the terms remaining AFTER pruning by\n"
"                          --parameters.\n"
"                        * You can use \"--tests all\" to include all terms.\n"
	       );
    help_print("linear\tdosage\tvif", &help_ctrl, 0,
"  --vif <max VIF>     : Set VIF threshold for --linear multicollinearity check\n"
"                        (default 50).\n"
	       );
    help_print("linear\tlogistic\txchr-model", &help_ctrl, 0,
"  --xchr-model <code> : Set the X chromosome --linear/--logistic model.\n"
"                        0 = skip sex and haploid chromosomes\n"
"                        1 (default) = add sex as a covariate on X chromosome\n"
"                        2 = code male genotypes 0/2 instead of 0/1\n"
"                        3 = test for interaction between genotype and sex\n"
	       );
    help_print("lasso\tlasso-select-covars", &help_ctrl, 0,
"  --lasso-select-covars [cov(s)...] : Subject some or all covariates to LASSO\n"
"                                      model selection.\n"
	       );
#ifndef STABLE_BUILD
    help_print("lasso\tlasso-lambda", &help_ctrl, 0,
"  --lasso-lambda <iters> [h2]       : Customize LASSO warm-start procedure.\n"
"                                      (h2 required if not used with --lasso.)\n"
	       );
#endif
    help_print("adjust\tgc\tlog10\tqq-plot", &help_ctrl, 0,
"  --adjust ['gc'] ['log10'] ['qq-plot'] : Report some multiple-testing\n"
"                                          corrections.\n"
	       );
    help_print("adjust\tlambda", &help_ctrl, 0,
"  --lambda <val>   : Set genomic control lambda for --adjust.\n"
	       );
    help_print("ci\tassoc\tlinear\tlogistic\tmh\tbd\ttdt", &help_ctrl, 0,
"  --ci <size>      : Report confidence intervals for odds ratios.\n"
	       );
    help_print("pfilter", &help_ctrl, 0,
"  --pfilter <val>  : Filter out association test results with higher p-values.\n"
	       );
    help_print("aperm", &help_ctrl, 0,
"  --aperm <min perms - 1> [max perms] [alpha] [beta] [init interval] [slope] :\n"
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
"      <slope>p + <init interval> more permutations (rounded down).  Default\n"
"      initial interval is 1, and default slope is 0.001.\n"
	       );
    help_print("mperm-save\tmperm-save-all", &help_ctrl, 0,
"  --mperm-save     : Save best max(T) permutation test statistics.\n"
"  --mperm-save-all : Save all max(T) permutation test statistics.\n"
	       );
    help_print("set-p\tset-r2\tset-max\tset-test-lambda\tset-test\twrite-set-r2\tset-r2-phase", &help_ctrl, 0,
"  --set-p <p-val>        : Adjust set test significant variant p-value ceiling\n"
"                           (default 0.05).\n"
"  --set-r2 [v] ['write'] : Adjust set test significant variant pairwise r^2\n"
"                           ceiling (default 0.5).  'write' causes violating\n"
"                          pairs to be dumped to <output prefix>.ldset.\n"
"  --set-max <ct>         : Adjust set test maximum # of significant variants\n"
"                           considered per set (default 5).\n"
"  --set-test-lambda <v>  : Specify genomic control correction for set test.\n"
	       );
    help_print("annotate\tborder\tannotate-snp-field", &help_ctrl, 0,
"  --border <kbs>            : Extend --annotate range intervals by given # kbs.\n"
"  --annotate-snp-field <nm> : Set --annotate variant ID field name.\n"
	       );
    help_print("clump-p1\tclump-p2\tclump-r2\tclump-kb\tclump-snp-field\tclump-field\tclump", &help_ctrl, 0,
"  --clump-p1 <pval> : Set --clump index var. p-value ceiling (default 1e-4).\n"
"  --clump-p2 <pval> : Set --clump secondary p-value threshold (default 0.01).\n"
"  --clump-r2 <r^2>  : Set --clump r^2 threshold (default 0.5).\n"
"  --clump-kb <kbs>  : Set --clump kb radius (default 250).\n"
"  --clump-snp-field <n...>  : Set --clump variant ID field name (default\n"
"                              'SNP').  With multiple field names, earlier names\n"
"                              take precedence over later ones.\n"
"  --clump-field <name...>   : Set --clump p-value field name (default 'P').\n"
	       );
    help_print("clump-allow-overlap\tclump", &help_ctrl, 0,
"  --clump-allow-overlap     : Let --clump non-index vars. join multiple clumps.\n"
	       );
    help_print("clump-verbose\tclump", &help_ctrl, 0,
"  --clump-verbose           : Request extended --clump report.\n"
	       );
    help_print("clump-annotate\tclump-verbose\tclump-best\tclump", &help_ctrl, 0,
"  --clump-annotate <hdr...> : Include named extra fields in --clump-verbose and\n"
"                              --clump-best reports.  (Field names can be\n"
"                              separated with spaces or commas.)\n"
	       );
    help_print("clump-range\tclump-range-border\tclump", &help_ctrl, 0,
"  --clump-range <filename>  : Report overlaps between clumps and regions.\n"
"  --clump-range-border <kb> : Stretch regions in --clump-range file.\n"
	       );
    help_print("clump-index-first\tclump-replicate\tclump", &help_ctrl, 0,
"  --clump-index-first       : Extract --clump index vars. from only first file.\n"
"  --clump-replicate         : Exclude clumps which contain secondary results\n"
"                              from only one file.\n"
	       );
    help_print("clump-best\tclump", &help_ctrl, 0,
"  --clump-best              : Report best proxy for each --clump index var.\n"
	       );
    help_print("meta-analysis-chr-field\tmeta-analysis-snp-field\tmeta-analysis-bp-field\tmeta-analysis-a1-field\tmeta-analysis-a2-field\tmeta-analysis-p-field\tmeta-analysis-se-field\tmeta-analysis-ess-field\tmeta-analysis", &help_ctrl, 0,
"  --meta-analysis-chr-field <n...> : Set --meta-analysis chromosome, variant\n"
"  --meta-analysis-snp-field <n...>   ID, position, A1/A2 allele, p-value,\n"
"  --meta-analysis-bp-field <n...>    standard error, and/or effective sample\n"
"  --meta-analysis-a1-field <n...>    size field names.\n"
"  --meta-analysis-a2-field <n...>    Defaults are 'CHR', 'SNP', 'BP', 'A1',\n"
"  --meta-analysis-p-field <n...>     'A2', 'P', 'SE', and 'NMISS',\n"
"  --meta-analysis-se-field <n...>    respectively.  When multiple parameters\n"
"  --meta-analysis-ess-field <n...>   are given to these flags, earlier names\n"
"                                     take precedence over later ones.\n"
"                                     Note that, if the numbers of cases and\n"
"                                     controls are unequal, effective sample\n"
"                                     size should be\n"
"                                       4 / (1/<# cases> + 1/<# controls>).\n"
	       );
    help_print("meta-analysis-report-dups\tmeta-analysis", &help_ctrl, 0,
"  --meta-analysis-report-dups      : When a variant appears multiple times in\n"
"                                     in the same file, report that.\n"
	       );
    help_print("gene-list-border\tgene-report\tgene-subset\tgene-list\tgene-report-snp-field", &help_ctrl, 0,
"  --gene-list-border <kbs>   : Extend --gene-report regions by given # of kbs.\n"
"  --gene-subset <filename>   : Specify gene name subset for --gene-report.\n"
"  --gene-report-snp-field <> : Set --gene-report variant ID field name (default\n"
"                               'SNP').  Only relevant with --extract.\n"
	       );
    help_print("fast-epistasis\tepistasis\tgap\tepi1\tepi2", &help_ctrl, 0,
"  --gap <kbs>      : Set \"--fast-epistasis case-only\" min. gap (default 1000).\n"
"  --epi1 <p-value> : Set --[fast-]epistasis reporting threshold (default\n"
"                     5e-6 for 'boost', 1e-4 otherwise).\n"
"  --epi2 <p-value> : Set threshold for contributing to SIG_E count (def. 0.01).\n"
	       );
    help_print("fast-epistasis\tje-cellmin", &help_ctrl, 0,
"  --je-cellmin <n> : Set required number of observations per 3x3x2 contingency\n"
"                     table cell for joint-effects test (default 5).\n"
	       );
    help_print("score\tq-score-file\tq-score-range", &help_ctrl, 0,
"  --q-score-range <range file> <data file> [i] [j] ['header'] :\n"
"    Apply --score to subset(s) of variants in the primary score list based\n"
"    on e.g. p-value ranges.\n"
"    * The first file should have range labels in the first column, p-value\n"
"      lower bounds in the second column, and upper bounds in the third column.\n"
"      Lines with too few entries, or nonnumeric values in the second or third\n"
"      column, are ignored.\n"
"    * The second file should contain a variant ID and a p-value on each\n"
"      nonempty line (except possibly the first).  Variant IDs are read from\n"
"      column #i and p-values are read from column #j, where i defaults to 1 and\n"
"      j defaults to i+1.  The 'header' modifier causes the first nonempty line\n"
"      of this file to be skipped.\n"
	       );
#if defined __cplusplus && !defined _WIN32
    help_print("R\tR-port\tR-host\tR-socket", &help_ctrl, 0,
"  --R-port <port #>  : Connect to Rserve on a port other than 6311.\n"
"  --R-host <host>    : Connect to Rserve host.\n"
"  --R-socket <sock>  : Connect to Rserve socket.\n"
	       );
#endif
    help_print("parallel\tgenome-lists", &help_ctrl, 0,
"  --parallel <k> <n> : Divide the output matrix into n pieces, and only compute\n"
"                       the kth piece.  The primary output file will have the\n"
"                       piece number included in its name, e.g. " PROG_NAME_STR ".rel.13 or\n"
"                       " PROG_NAME_STR ".rel.13.gz if k is 13.  Concatenating these files\n"
"                       in order will yield the full matrix of interest.  (Yes,\n"
"                       this can be done before unzipping.)\n"
"                       N.B. This generally cannot be used to directly write a\n"
"                       symmetric square matrix.  Choose square0 or triangle\n"
"                       shape instead, and postprocess as necessary.\n"
	       );
    help_print("memory", &help_ctrl, 0,
"  --memory <val>     : Set size, in MB, of initial workspace malloc attempt.\n"
"                       (Practically mandatory when using GNU parallel.)\n"
	       );
    help_print("threads\tthread-num\tnum_threads", &help_ctrl, 0,
"  --threads <val>    : Set maximum number of concurrent threads.\n"
"                       This has one known limitation: some BLAS/LAPACK linear\n"
"                       algebra operations are multithreaded in a way that PLINK\n"
"                       cannot control.  If this is problematic, you should\n"
"                       recompile against single-threaded BLAS/LAPACK.\n"
	       );
    help_print("d\tsnps", &help_ctrl, 0,
"  --d <char>         : Change variant/covariate range delimiter (normally '-').\n"
	       );
    help_print("seed", &help_ctrl, 0,
"  --seed <val...>    : Set random number seed(s).  Each value must be an\n"
"                       integer between 0 and 4294967295 inclusive.\n"
	       );
    help_print("perm-batch-size", &help_ctrl, 0,
"  --perm-batch-size <val> : Set number of permutations per batch for some\n"
"                            permutation tests.\n"
	       );
    help_print("output-min-p", &help_ctrl, 0,
"  --output-min-p <p> : Specify minimum p-value to write to reports.\n"
	       );
    help_print("debug", &help_ctrl, 0,
"  --debug            : Use slower, more crash-resistant logging method.\n"
	       );
#ifndef STABLE_BUILD
    if (!param_ct) {
      fputs(
"\nThese flags only apply to .cnv fileset analysis:\n"
, stdout);
    }
    help_print("cnv-del\tcnv-dup", &help_ctrl, 0,
"  --cnv-del                : Exclude all variants with more than one copy.\n"
"  --cnv-dup                : Exclude all variants with fewer than three copies.\n"
	       );
    help_print("cnv-kb\tcnv-max-kb", &help_ctrl, 0,
"  --cnv-kb <kb len>        : Exclude segments shorter than the given length.\n"
"  --cnv-max-kb <kb len>    : Exclude segments longer than the given length.\n"
	       );
    help_print("cnv-score\tcnv-max-score", &help_ctrl, 0,
"  --cnv-score <val>        : Exclude all variants with confidence score < val.\n"
"  --cnv-max-score <val>    : Exclude all variants with confidence score > val.\n"
	       );
    help_print("cnv-sites\tcnv-max-sites", &help_ctrl, 0,
"  --cnv-sites <ct>         : Exclude all segments with fewer than ct probes.\n"
"  --cnv-max-sites <ct>     : Exclude all segments with more than ct probes.\n"
	       );
    help_print("cnv-intersect\tcnv-exclude\tcnv-subset\tcnv-overlap\tcnv-region-overlap\tcnv-union-overlap\tcnv-disrupt", &help_ctrl, 0,
"  --cnv-intersect <fname>  : Only include segments which intersect a region in\n"
"                             the given region list.\n"
"  --cnv-exclude <fname>    : Exclude all segments which intersect a region in\n"
"                             the given region list.\n"
"  --cnv-subset <fname>     : Ignore all regions in the --cnv-intersect/-exclude\n"
"                             /-count list that aren't named in the given file.\n"
"  --cnv-overlap <x>        : Only count intersections of length at least xn,\n"
"                             where n is the segment size.\n"
"  --cnv-region-overlap <x> : x >= <overlap> / <region size>.\n"
"  --cnv-union-overlap <x>  : x >= <overlap> / <union size>.\n"
"  --cnv-disrupt            : Only include/exclude segments with an endpoint in\n"
"                             a region.\n"
	       );
    help_print("cnv-freq-exclude-above\tcnv-freq-exclude-below\tcnv-freq-exclude-exact\tcnv-freq-include-exact\tcnv-freq-overlap\tcnv-freq-method2", &help_ctrl, 0,
"  --cnv-freq-exclude-above <k> : Exclude all segments where any portion is\n"
"                                 included by more than k total segments.\n"
"  --cnv-freq-exclude-below <k> : Exclude all segments where no portion is\n"
"                                 included by k or more total segments.\n"
"  --cnv-freq-exclude-exact <k> : Exclude all segments which have a portion\n"
"                                 included by at least k total segments, but no\n"
"                                 portion included by more.\n"
"  --cnv-freq-include-exact <k> : Reverse of --cnv-freq-exclude-exact.\n"
"  --cnv-freq-overlap [x]   : Only count portions of length at least xn, where n\n"
"                             is the segment size.\n"
"  --cnv-freq-method2 [x]   : Causes k to instead be compared against the number\n"
"                             of segments for which x >= <overlap> / <union>.\n"
	       );
    help_print("cnv-exclude-off-by-1", &help_ctrl, 0,
"  --cnv-exclude-off-by-1   : Exclude .cnv segments where the terminal .cnv.map\n"
"                             entry is off by 1.\n"
	       );
    help_print("cnv-test-window\tcnv-test", &help_ctrl, 0,
"  --cnv-test-window <size> : Specify window size (in kb) for CNV assoc. test.\n"
	       );
    help_print("cnv-count\tcnv-indiv-perm\tcnv-enrichment-test", &help_ctrl, 0,
"  --cnv-count <fname>      : Specify region list for --cnv-indiv-perm\n"
"                             (optional) or --cnv-enrichment-test (required).\n"
	       );
#endif
    if (!param_ct) {
      fputs(
"\nPrimary methods paper:\n"
"Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ (2015)\n"
"Second-generation PLINK: rising to the challenge of larger and richer datasets.\n"
"GigaScience, 4.\n\n"
"For further documentation and support, consult the main webpage\n"
"(https://www.cog-genomics.org/plink/1.9 ) and/or the mailing list\n"
"(https://groups.google.com/d/forum/plink2-users ).\n"
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
	  if (help_ctrl.param_lens[arg_uidx] + col_num > 76) {
	    putc_unlocked('\n', stdout);
	    col_num = 2 + help_ctrl.param_lens[arg_uidx];
	  } else {
	    putc_unlocked(' ', stdout);
	    col_num += 3 + help_ctrl.param_lens[arg_uidx];
	  }
	  putc_unlocked('\'', stdout);
	  fputs(argv[arg_uidx], stdout);
	  putc_unlocked('\'', stdout);
	} else {
	  if (help_ctrl.param_lens[arg_uidx] + col_num > 75) {
	    putc_unlocked('\n', stdout);
	    col_num = 3 + help_ctrl.param_lens[arg_uidx];
	  } else {
	    putc_unlocked(' ', stdout);
	    col_num += 4 + help_ctrl.param_lens[arg_uidx];
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
	putc_unlocked((help_ctrl.param_lens[arg_uidx] + col_num > 75)? '\n' : ' ', stdout);
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
      retval = RET_NOMEM;
    }
    free_cond(help_ctrl.param_lens);
    free_cond(help_ctrl.all_match_arr);
    if (help_ctrl.argv && (help_ctrl.argv != argv)) {
      free(help_ctrl.argv);
    }
  }
  return retval;
}

#include "wdist_common.h"

int32_t edit1_match(int32_t len1, char* s1, int32_t len2, char* s2) {
  // permit one difference of the following forms:
  // - inserted/deleted character
  // - replaced character
  // - adjacent pair of swapped characters
  int32_t diff_found = 0;
  int32_t pos = 0;
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

#define MAX_EQUAL_HELP_PARAMS 16

typedef struct {
  int32_t iters_left;
  uint32_t param_ct;
  char** argv;
  uintptr_t unmatched_ct;
  uintptr_t* all_match_arr;
  uintptr_t* prefix_match_arr;
  uintptr_t* perfect_match_arr;
  uint32_t* param_lens;
  int32_t preprint_newline;
} Help_ctrl;

void help_print(const char* cur_params, Help_ctrl* help_ctrl_ptr, int32_t postprint_newline, const char* payload) {
  // unmatched_ct fixed during call, *unmatched_ct_ptr may decrease
  uint32_t unmatched_ct = help_ctrl_ptr->unmatched_ct;
  int32_t print_this = 0;
  int32_t cur_param_lens[MAX_EQUAL_HELP_PARAMS];
  char* cur_param_start[MAX_EQUAL_HELP_PARAMS];
  int32_t cur_param_ct;
  int32_t cur_param_idx;
  uint32_t arg_uidx;
  uint32_t arg_idx;
  int32_t uii;
  int32_t payload_len;
  char* payload_ptr;
  char* line_end;
  char* payload_end;
  if (help_ctrl_ptr->param_ct) {
    strcpy(tbuf, cur_params);
    cur_param_ct = 1;
    cur_param_start[0] = tbuf;
    payload_ptr = strchr(tbuf, '\t');
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
	    arg_uidx = next_non_set_unsafe(help_ctrl_ptr->all_match_arr, arg_uidx);
	    for (cur_param_idx = 0; cur_param_idx < cur_param_ct; cur_param_idx++) {
	      if (!strcmp(cur_param_start[cur_param_idx], help_ctrl_ptr->argv[arg_uidx])) {
		set_bit_noct(help_ctrl_ptr->perfect_match_arr, arg_uidx);
		set_bit_noct(help_ctrl_ptr->prefix_match_arr, arg_uidx);
		set_bit_sub(help_ctrl_ptr->all_match_arr, arg_uidx, &(help_ctrl_ptr->unmatched_ct));
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
	    arg_uidx = next_non_set_unsafe(help_ctrl_ptr->all_match_arr, arg_uidx);
	    uii = help_ctrl_ptr->param_lens[arg_uidx];
	    for (cur_param_idx = 0; cur_param_idx < cur_param_ct; cur_param_idx++) {
	      if (cur_param_lens[cur_param_idx] > uii) {
		if (!memcmp(help_ctrl_ptr->argv[arg_uidx], cur_param_start[cur_param_idx], uii)) {
		  set_bit_noct(help_ctrl_ptr->prefix_match_arr, arg_uidx);
		  set_bit_sub(help_ctrl_ptr->all_match_arr, arg_uidx, &(help_ctrl_ptr->unmatched_ct));
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
	if (is_set(help_ctrl_ptr->prefix_match_arr, arg_uidx)) {
	  if (!print_this) {
	    if (is_set(help_ctrl_ptr->perfect_match_arr, arg_uidx)) {
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
	      set_bit_sub(help_ctrl_ptr->all_match_arr, arg_uidx, &(help_ctrl_ptr->unmatched_ct));
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
	  putchar('\n');
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
	  memcpyx(tbuf, payload_ptr, uii, 0);
	  fputs(tbuf, stdout);
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
  uint32_t param_ctl = (param_ct + (BITCT - 1)) / BITCT;
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
  help_ctrl.param_lens = NULL;
  help_ctrl.all_match_arr = NULL;
  help_ctrl.argv = NULL;
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
    fill_ulong_zero(help_ctrl.all_match_arr, param_ctl * 3);
    help_ctrl.prefix_match_arr = &(help_ctrl.all_match_arr[param_ctl]);
    help_ctrl.perfect_match_arr = &(help_ctrl.all_match_arr[param_ctl * 2]);
    help_ctrl.preprint_newline = 1;
  } else {
    help_ctrl.argv = NULL;
    fputs(
"\nIn the command line flag definitions that follow,\n"
"  * [square brackets] denote a required parameter, where the text between the\n"
"    brackets describes its nature.\n"
"  * <angle brackets> denote an optional modifier (or if '|' is present, a set\n"
"    of mutually exclusive optional modifiers).  Use the EXACT text in the\n"
"    definition, e.g. '--dummy acgt'.\n"
"  * There's one exception to the angle brackets/exact text rule: when an angle\n"
"    bracket term ends with '=[value]', '[value]' designates a variable\n"
"    parameter.\n"
"  * {curly braces} denote an optional parameter, where the text between the\n"
"    braces describes its nature.\n"
"  * An ellipsis (...) indicates that you may enter multiple parameters of the\n"
"    specified type.\n"
, stdout);
    fputs(cmdline_format_str, stdout);
    fputs(
"Each WDIST run requires exactly one main input fileset.  The following flags\n"
"are available for defining its form and location:\n\n"
, stdout);
  }
  do {
    help_print("bfile\tbed\tbim\tfam", &help_ctrl, 1,
"  --bfile {prefix} : Specify .bed + .bim + .fam prefix (default 'wdist').\n"
"  --bed [filename] : Specify full name of .bed file.\n"
"  --bim [filename] : Specify full name of .bim file.\n"
"  --fam [filename] : Specify full name of .fam file.\n\n"
	       );
    help_print("file\tped\tmap", &help_ctrl, 1,
"  --file {prefix}  : Specify .ped + .map filename prefix (default 'wdist').\n"
"  --ped [filename] : Specify full name of .ped file.\n"
"  --map [filename] : Specify full name of .map file.\n\n"
	       );
    help_print("tfile\ttped\ttfam", &help_ctrl, 1,
"  --tfile {prefix} : Specify .tped + .tfam filename prefix (default 'wdist').\n"
"  --tped [fname]   : Specify full name of .tped file.\n"
"  --tfam [fname]   : Specify full name of .tfam file.\n\n"
	       );
    help_print("lfile", &help_ctrl, 1,
"  --lfile {prefix} : Specify .lgen + .map + .fam (long-format fileset) prefix.\n\n"
	       );
    help_print("data\tgen\tsample", &help_ctrl, 1,
"  --data {prefix}  : Specify Oxford .gen + .sample prefix (default 'wdist').\n"
"  --gen [filename] : Specify full name of .gen file.\n"
"  --sample [fname] : Specify full name of .sample file.\n\n"
	       );
    help_print("cfile\tcnv-list\tgfile", &help_ctrl, 1,
"  --cfile [prefix] : Specify .cnv + .fam + .cnv.map (segmental CNV) prefix.\n"
"  --cnv-list [fn]  : Specify full name of .cnv file.\n"
"  --gfile [prefix] : Specify .gvar + .fam + .map (genetic variant) prefix.\n\n"
	       );
    help_print("grm\tgrm-bin\trel-cutoff\tgrm-cutoff", &help_ctrl, 1,
"  --grm {prefix}   : Specify .grm.gz + .grm.id (GCTA rel. matrix) prefix.\n"
"  --grm-bin {prfx} : Specify .grm.bin + .grm.N.bin + .grm.id (GCTA triangular\n"
"                     binary relationship matrix) filename prefix.\n\n"
	       );
    help_print("dummy", &help_ctrl, 1,
"  --dummy [indiv ct] [marker ct] {missing geno freq} {missing pheno freq}\n"
"          <acgt | 1234 | 12> <scalar-pheno>\n"
"    This generates a fake input dataset with the specified number of\n"
"    individuals and markers.  By default, the missing genotype and phenotype\n"
"    frequencies are zero, and genotypes are As and Bs (change the latter with\n"
"    'acgt'/'1234'/'12').  The 'scalar-pheno' modifier causes a normally\n"
"    distributed scalar phenotype to be generated instead of a binary one.\n\n"
	       );
    help_print("simulate\tsimulate-qt", &help_ctrl, 1,
"  --simulate [simulation parameter file] <tags | haps> <acgt | 1234 | 12>\n"
"  --simulate-qt [simulation parameter file] <tags | haps> <acgt | 1234 | 12>\n"
"    --simulate generates a fake input dataset with disease-associated SNPs,\n"
"    while --simulate-qt generates a dataset with quantitative trait loci.\n\n"
	       );
    if (!param_ct) {
      fputs(
"Output files have names of the form 'wdist.{extension}' by default.  You can\n"
"change the 'wdist' prefix with\n\n"
, stdout);
    }
    help_print("out", &help_ctrl, 1,
"  --out [prefix]   : Specify prefix for output files.\n\n"
	       );
    if (!param_ct) {
      fputs(
"WDIST automatically converts PLINK text filesets to binary during the loading\n"
"process (the new fileset is saved to {output prefix}.bed + .bim + .fam, unless\n"
"that would conflict with the output of another command like --make-bed).  You\n"
"are encouraged to directly use the new binary fileset in future runs.\n\n"
"Every run also requires at least one of the following commands (unless you just\n"
"want automatic text-to-binary conversion or the default CNV report):\n\n"
, stdout);
    }
    help_print("make-bed", &help_ctrl, 1,
"  --make-bed\n"
"    Creates a new binary fileset.  Unlike the automatic text-to-binary\n"
"    converter (which only respects --autosome and --chr), this supports all of\n"
"    WDIST's filtering flags.\n"
	       );
    help_print("recode\trecode12\ttab\ttranspose\trecode-lgen\trecodeAD\trecodead\trecodeA\trecodea\trecode-rlist\trecode-allele\tlist\twith-reference\trecode-vcf\tfid\tiid", &help_ctrl, 1,
"  --recode <12> <compound-genotypes> <23 | A | AD | lgen | lgen-ref | list |\n"
"           rlist | transpose | vcf | vcf-fid | vcf-iid> <tab | tabx | spacex>\n"
"    Creates a new text fileset with all filters applied.\n"
"    * The '12' modifier causes all alleles to be coded as 1s and 2s.\n"
"    * The 'compound-genotypes' modifier removes the space between pairs of\n"
"      genotype codes for the same marker.\n"
"    * The '23' modifier causes a 23andMe-formatted file to be generated.  This\n"
"      can only be used on a single individual's data (--keep may be handy).\n"
"    * The 'AD' modifier causes an additive + dominant component file, suitable\n"
"      for loading from R, to be generated instead.  If you don't want the\n"
"      dominant component, use 'A' instead.\n"
"    * The 'lgen' modifier causes a long-format fileset to be generated instead,\n"
"      (loadable with --lfile), while 'lgen-ref' generates a (usually) smaller\n"
"      long-format fileset loadable with --lfile + --reference.\n"
"    * The 'list' modifier creates a genotype-based list, while 'rlist'\n"
"      creates a rare-genotype fileset.  (These formats are not directly\n"
"      loadable by WDIST or PLINK.)\n"
"    * 'transpose' creates a transposed text fileset (loadable with --tfile).\n"
"    * 'vcf', 'vcf-fid', and 'vcf-iid' result in production of a VCFv4.0 file.\n"
"      'vcf-fid' and 'vcf-iid' cause family IDs or within-family IDs\n"
"      respectively to be used for the sample IDs in the last header row, while\n"
"      'vcf' merges both IDs and puts an underscore between them.\n"
"    * The 'tab' modifier makes the output mostly tab-delimited instead of\n"
"      mostly space-delimited.  'tabx' and 'spacex' force all tabs and all\n"
"      spaces, respectively.\n\n"
	       );
    help_print("write-covar", &help_ctrl, 1,
"  --write-covar\n"
"    If a --covar file is loaded, this creates a revised covariate file after\n"
"    applying all filters.  (This automatically happens if --make-bed or\n"
"    --recode is specified.)\n\n"
	       );
    help_print("merge\tbmerge\tmerge-list\tmerge-mode", &help_ctrl, 1,
"  --merge [.ped filename] [.map filename]\n"
"  --merge [text fileset prefix]\n"
"  --bmerge [.bed filename] [.bim filename] [.fam filename]\n"
"  --bmerge [binary fileset prefix]\n"
"    Merges the given fileset with the initially loaded fileset.  If you specify\n"
"    --make-bed or '--recode lgen', the initial merge result is written to\n"
"    {output prefix}-merge.bed + .bim + .fam, filtering is performed, and then\n"
"    the post-filtering data is written to {output prefix}.bed + .bim + .fam.\n"
"    Otherwise, the merged data is written directly to\n"
"    {output prefix}.bed + .bim + .fam.\n"
"  --merge-list [filename]\n"
"    Merge all filesets named in the text file with the initially loaded\n"
"    fileset.  The text file is interpreted as follows:\n"
"    * If a line contains only one name, it is assumed to be the prefix for a\n"
"      binary fileset.\n"
"    * If a line contains exactly two names, they are assumed to be the full\n"
"      filenames for a text fileset (.ped first, then .map).\n"
"    * If a line contains exactly three names, they are assumed to be the full\n"
"      filenames for a binary fileset (.bed, then .bim, then .fam).\n\n"
	       );
    help_print("write-snplist", &help_ctrl, 1,
"  --write-snplist\n"
"    Writes a .snplist file listing the names of all markers that pass the\n"
"    filters and inclusion thresholds you've specified.\n\n"
	       );
    help_print("freq\tfreqx\tfrqx\tcounts", &help_ctrl, 1,
"  --freq <counts>\n"
"  --freqx\n"
"    --freq generates an allele frequency report identical to that of PLINK\n"
"    --freq (or --freq --counts, if the 'counts' modifier is used).  Using the\n"
"    --freqx flag instead causes the MAF and NCHROBS columns to be replaced with\n"
"    homozygote and heterozygote counts, which (when reloaded with --read-freq)\n"
"    allow distance matrix terms to be weighted consistently through multiple\n"
"    filtering runs.\n\n"
		);
    help_print("hardy", &help_ctrl, 1,
"  --hardy\n"
"    Generates a Hardy-Weinberg exact test p-value report.  (This does NOT\n"
"    filter on the p-value; use --hwe for that.)\n\n"
	       );
    help_print("ibc\thet", &help_ctrl, 1,
"  --ibc\n"
"    Calculates inbreeding coefficients in three different ways.\n"
"    * For more details, see Yang J, Lee SH, Goddard ME and Visscher PM.  GCTA:\n"
"      a tool for Genome-wide Complex Trait Analysis.  Am J Hum Genet. 2011 Jan\n"
"      88(1): 76-82.  This paper also describes the relationship matrix\n"
"      computation we implement.\n\n"
	       );
    help_print("distance\tregress-pcs-distance", &help_ctrl, 1,
"  --distance <square | square0 | triangle> <gz | bin> <ibs> <1-ibs> <allele-ct>\n"
"             <3d> <flat-missing>\n"
"    Writes a lower-triangular tab-delimited table of (weighted) genomic\n"
"    distances in allele count units to {output prefix}.dist, and a list of the\n"
"    corresponding family/individual IDs to {output prefix}.dist.id.  The first\n"
"    row of the .dist file contains a single {genotype 1-genotype 2} distance,\n"
"    the second row has the {genotype 1-genotype 3} and {genotype 2-genotype 3}\n"
"    distances in that order, etc.\n"
"    * If the 'square' or 'square0' modifier is present, a square matrix is\n"
"      written instead; 'square0' fills the upper right triangle with zeroes.\n"
"    * If the 'gz' modifier is present, a compressed .dist.gz file is written\n"
"      instead of a plain text file.\n"
"    * If the 'bin' modifier is present, a binary (square) matrix of\n"
"      double-precision floating point values, suitable for loading from R, is\n"
"      instead written to {output prefix}.dist.bin.  This can be combined with\n"
"      'square0' if you still want the upper right zeroed out, or 'triangle' if\n"
"      you don't want to pad the upper right at all.\n"
"    * If the 'ibs' modifier is present, an identity-by-state matrix is written\n"
"      to {output prefix}.mibs.  '1-ibs' causes distances expressed as genomic\n"
"      proportions (i.e. 1 - IBS) to be written to {output prefix}.mdist.\n"
"      Combine with 'allele-ct' if you want to generate the usual .dist file as\n"
"      well.\n"
"    * With dosage data, the '3d' modifier considers 0-1-2 allele count\n"
"      probabilities separately, instead of collapsing them into an expected\n"
"      value and a missingness probability.\n"
"    * By default, distance rescaling in the presence of missing markers is\n"
"      sensitive to allele count distributions: if allele A contributes, on\n"
"      average, twice as much to other pairwise distances as allele B, a missing\n"
"      allele A will result in twice as large of a missingness correction.  To\n"
"      turn this off, use the 'flat-missing' modifier.\n"
	       );
    help_print("matrix\tdistance-matrix", &help_ctrl, 1,
"  --matrix\n"
"  --distance-matrix\n"
"    These generate space-delimited text matrices, and are included for\n"
"    backwards compatibility with old scripts relying on the corresponding PLINK\n"
"    flags.  New scripts should migrate to '--distance ibs' and '--distance\n"
"    1-ibs', which support output formats better suited to parallel computation\n"
"    and have more accurate handling of missing markers.\n\n"
		);
    help_print("genome\tZ-genome\tgenome-full\tunbounded", &help_ctrl, 1,
"  --genome <gz> <full> <unbounded>\n"
"    Identity-by-descent analysis.  This yields the same output as PLINK\n"
"    --genome/--Z-genome, and the 'full' and 'unbounded' modifiers have the same\n"
"    effect as PLINK's --genome-full and --unbounded flags.\n\n"
		);
    help_print("assoc\tmodel\tfisher\tperm\tmperm\tperm-count\tcounts\tp2\tmodel-dom\tmodel-gen\tmodel-rec\tmodel-trend\tgenedrop\tqt-means\ttrend", &help_ctrl, 1,
"  --assoc <perm | mperm=[value]> <genedrop> <perm-count> <fisher> <counts> <p2>\n"
"  --assoc <perm | mperm=[value]> <perm-count> <qt-means> <lin>\n"
"  --model <perm | mperm=[value]> <genedrop> <perm-count> <fisher | trend-only>\n"
"          <dom | rec | gen | trend>\n"
"    Basic association analysis report.\n"
"    Given a case/control phenotype, --assoc performs a 1df chi-square allelic\n"
"    test, while --model performs 4 other tests as well (1df dominant gene\n"
"    action, 1df recessive gene action, 2df genotypic, Cochran-Armitage trend).\n"
"    * With 'fisher', Fisher's exact test is used to generate p-values.\n"
"    * 'perm' causes an adaptive permutation test to be performed.\n"
"    * 'mperm=[value]' causes a max(T) permutation test with the specified\n"
"      number of replications to be performed.\n"
"    * 'genedrop' causes offspring genotypes to be regenerated via gene-dropping\n"
"      in the permutation test.\n"
"    * 'perm-count' causes the permutation test report to include counts instead\n"
"      of frequencies.\n"
"    * 'counts' causes --assoc to report allele counts instead of frequencies.\n"
"    * 'p2' changes the --assoc permutation test used (see PLINK documentation).\n"
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
    help_print("gxe\tmcovar", &help_ctrl, 1,
"  --gxe {covariate index}\n"
"    Given both a quantitative phenotype and a dichotomous covariate loaded with\n"
"    --covar defining two groups, --gxe compares the regression coefficient\n"
"    derived from considering only members of one group to the regression\n"
"    coefficient derived from considering only members of the other.  By\n"
"    default, the first covariate in the --covar file defines the groups; use\n"
"    e.g. '--gxe 3' to base them on the third covariate instead.\n\n"
	       );
    help_print(
"indep\tindep-pairwise", &help_ctrl, 1,
"  --indep [window size]<kb> [step size (markers)] [VIF threshold]\n"
"  --indep-pairwise [window size]<kb> [step size (markers)] [r^2 threshold]\n"
"    Generates a list of markers in approximate linkage equilibrium.  With the\n"
"    'kb' modifier, the window size is in kilobase units instead of marker\n"
"    count.  (Pre-'kb' space is optional, i.e. '--indep-pairwise 500 kb 5 0.5'\n"
"    and '--indep-pairwise 500kb 5 0.5' have the same effect.)\n"
"    Note that you need to rerun WDIST using --extract or --exclude on the\n"
"    .prune.in/.prune.out file to apply the list to another computation.\n\n"
		);
    help_print("make-rel", &help_ctrl, 1,
"  --make-rel <square | square0 | triangle> <gz | bin> <cov | ibc2 | ibc3>\n"
"             <single-prec>\n"
"    Writes a lower-triangular variance-standardized relationship (coancestry)\n"
"    matrix to {output prefix}.rel, and corresponding IDs to\n"
"    {output prefix}.rel.id.\n"
"    * 'square', 'square0', 'triangle', 'gz', and 'bin' act as they do on\n"
"      --distance.\n"
"    * The 'cov' modifier removes the variance standardization step, causing a\n"
"      covariance matrix to be calculated instead.\n"
"    * By default, the diagonal elements in the relationship matrix are based on\n"
"      --ibc's Fhat1; use the 'ibc2' or 'ibc3' modifiers to base them on Fhat2\n"
"      or Fhat3 instead.\n"
"    * WDIST normally performs this calculation with double-precision floating\n"
"      point numbers.  The 'single-prec' modifier switches to single-precision\n"
"      arithmetic, which is generally good enough; this decreases memory usage\n"
"      and speeds up computation.\n"
               );
    help_print("make-grm\tmake-grm-bin\tgrm\tgrm-bin", &help_ctrl, 1,
"  --make-grm <no-gz> <cov | ibc2 | ibc3> <single-prec>\n"
"  --make-grm-bin <cov | ibc2 | ibc3>\n"
"    --make-grm writes the relationships in GCTA's gzipped list format, which\n"
"    describes one pair per line, while --make-grm-bin writes them in GCTA's\n"
"    single-precision triangular binary format.  Note that these formats\n"
"    explicitly store the number of valid observations (where neither individual\n"
"    has a missing call) for each pair, which is useful input for some scripts.\n\n"
	       );
    help_print("rel-cutoff\tgrm-cutoff", &help_ctrl, 1,
"  --rel-cutoff {val}\n"
"  --grm-cutoff {val}\n"
"    Excludes one member of each pair of individuals with relatedness greater\n"
"    than the given cutoff value (default 0.025).  If no later operation will\n"
"    cause the list of remaining individuals to be written to disk, this will\n"
"    save it to {output prefix}.rel.id.\n"
"    Note that maximizing the remaining sample size is equivalent to the NP-hard\n"
"    maximum independent set problem, so we use a greedy algorithm instead of\n"
"    guaranteeing optimality.  (Use the --make-rel and --keep/--remove flags if\n"
"    you want to try to do better.)\n\n"
	       );
    help_print("regress-distance", &help_ctrl, 1,
"  --regress-distance {iters} {d}\n"
"    Linear regression of pairwise genomic distances on pairwise average\n"
"    phenotypes and vice versa, using delete-d jackknife for standard errors.\n"
"    Scalar phenotype data is required.\n"
"    * With less than two parameters, d is set to {number of people}^0.6 rounded\n"
"      down.  With no parameters, 100k iterations are run.\n\n"
	       );
    help_print("regress-pcs\tdistance\tregress-pcs-distance", &help_ctrl, 1,
"  --regress-pcs [.evec or .eigenvec filename] <normalize-pheno> <sex-specific>\n"
"                <clip> {max PCs}\n"
"    Linear regression of phenotypes and genotypes on the given list of\n"
"    principal components (produced by SMARTPCA or GCTA).  Output is a .gen +\n"
"    .sample fileset in the Oxford IMPUTE/SNPTEST v2 format.\n"
"    * The 'normalize-pheno' modifier converts all phenotype residuals to\n"
"      Z-scores.  When combined with 'sex-specific', the Z-scores are evaluated\n"
"      separately by sex.\n"
"    * The 'clip' modifier clips out-of-range genotype residuals.  Without it,\n"
"      they are represented as negative probabilities in the .gen file, which\n"
"      are invalid input for some programs.\n"
"    * By default, principal components beyond the 20th are ignored; change this\n"
"      by setting the max PCs parameter.\n"
"  --regress-pcs-distance [.evec/.eigenvec file] <normalize-pheno>\n"
"                         <sex-specific> {max PCs} <square | square0 | triangle>\n"
"                         <gz | bin> <ibs> <1-ibs> <allele-ct> <3d>\n"
"                         <flat-missing>\n"
"    High-speed combination of --regress-pcs and --distance (no .gen + .sample\n"
"    fileset is written to disk).\n\n"
	       );
    help_print("regress-rel", &help_ctrl, 1,
"  --regress-rel {iters} {d}\n"
"    Linear regression of pairwise genomic relationships on pairwise average\n"
"    phenotypes, and vice versa.  Defaults for iters and d are the same as for\n"
"    --regress-distance.\n\n"
	       );
    help_print("ibs-test\tgroupdist", &help_ctrl, 1,
"  --ibs-test {permutation count}\n"
"  --groupdist {iters} {d}\n"
"    Given dichotomous phenotype data, these commands consider three subsets of\n"
"    the distance matrix: pairs of affected individuals, affected-unaffected\n"
"    pairs, and pairs of unaffected individuals.  Each of these subsets has a\n"
"    distribution of pairwise genomic distances; --ibs-test uses permutation to\n"
"    estimate p-values re: which types of pairs are most similar, while\n"
"    --groupdist focuses on the differences between the centers of these\n"
"    distributions and estimates standard errors via delete-d jackknife.\n\n"
	       );
#ifndef NOLAPACK
    help_print("unrelated-heritability", &help_ctrl, 1,
"  --unrelated-heritability <strict> {tol} {initial covg} {initial covr}\n"
"    REML estimate of additive heritability, iterating with an accelerated\n"
"    variant of the EM algorithm until the rate of change of the log likelihood\n"
"    function is less than tol.  Scalar phenotype data is required.\n"
"    * The 'strict' modifier forces regular EM to be used.  tol defaults to\n"
"      10^{-7}, genomic covariance prior defaults to 0.45, and residual\n"
"      covariance prior defaults to (1 - covg).\n"
"    * For more details, see Vattikuti S, Guo J, Chow CC (2012) Heritability and\n"
"      Genetic Correlations Explained by Common SNPs for Metabolic Syndrome\n"
"      Traits.  PLoS Genet 8(3): e1002637.  doi:10.1371/journal.pgen.1002637\n\n"
	       );
#endif
    help_print("cnv-make-map", &help_ctrl, 1,
"  --cnv-make-map <short>\n"
"    Given a .cnv file, this generates the corresponding .cnv.map file needed\n"
"    by WDIST and PLINK's other CNV analysis commands.  The 'short' modifier\n"
"    causes entries needed by PLINK but not WDIST to be omitted.  (Now\n"
"    automatically invoked, with 'short', when necessary.)\n\n"
	       );
    help_print("cnv-check-no-overlap", &help_ctrl, 1,
"  --cnv-check-no-overlap\n"
"    Given a .cnv fileset, this checks for within-individual CNV overlaps.\n\n"
	       );
    help_print("cnv-write", &help_ctrl, 1,
"  --cnv-write <freq>\n"
"    Writes a new .cnv fileset, after applying all requested filters.  The\n"
"    'freq' modifier (which must be used with --cnv-freq-method2) causes an\n"
"    additional \"FREQ\" field to be written with CNV-CNV overlap counts.\n\n"
	       );
    help_print("cnv-indiv-perm\tcnv-test\tcnv-test-region\tcnv-enrichment-test\tmperm\tcnv-test-1sided\tcnv-test-2sided", &help_ctrl, 1,
"  --cnv-indiv-perm [permutation count]\n"
"  --cnv-test <1sided | 2sided> [permutation count]\n"
"  --cnv-test-region [permutation count]\n"
"  --cnv-enrichment-test {permutation count}\n"
"    Given a .cnv fileset,\n"
"    * --cnv-indiv-perm performs a case/control CNV burden test.\n"
"    * --cnv-test performs a basic permutation-based CNV association test.  By\n"
"      default, this is 1-sided for case/control phenotypes and 2-sided for\n"
"      quantitative traits; you can use '1sided'/'2sided' to force the opposite.\n"
"    * --cnv-test-region reports case/ctrl association test results by region.\n"
"    * --cnv-enrichment-test performs Raychaudhuri et al.'s geneset enrichment\n"
"      test.  Gene locations must be loaded with --cnv-count.\n\n"
	       );
    if (!param_ct) {
      fputs(
"The following other flags are supported.  (Order of operations is described at\n"
"https://www.cog-genomics.org/wdist/order .)\n"
, stdout);
    }
    help_print("script", &help_ctrl, 0,
"  --script [fname] : Include command-line options from file.\n"
	       );
    help_print("rerun", &help_ctrl, 0,
"  --rerun {log}    : Rerun commands in log (default 'wdist.log').\n"
	       );
    help_print("no-fid", &help_ctrl, 0,
"  --no-fid         : .fam/.ped file does not contain column 1 (family ID).\n"
	       );
    help_print("no-parents", &help_ctrl, 0,
"  --no-parents     : .fam/.ped file does not contain columns 3-4 (parents).\n"
	       );
    help_print("no-sex", &help_ctrl, 0,
"  --no-sex         : .fam/.ped file does not contain column 5 (sex).\n"
	       );
    help_print("no-pheno", &help_ctrl, 0,
"  --no-pheno       : .fam/.ped file does not contain column 6 (phenotype).\n"
	       );
    help_print("load-dists\tgroupdist\tregress-distance", &help_ctrl, 0,
"  --load-dists [f] : Load a binary TRIANGULAR distance matrix for --groupdist\n"
"                     or --regress-distance analysis, instead of recalculating\n"
"                     it from scratch.\n"
	       );
    help_print("silent", &help_ctrl, 0,
"  --silent         : Suppress output to console.\n"
	       );
    help_print("merge\tbmerge\tmerge-list\tmerge-mode", &help_ctrl, 0,
"  --merge-mode [n] : Adjust --merge/--bmerge/--merge-list behavior based on a\n"
"                     numeric code.\n"
"                     1 (default) = difference -> missing\n"
"                     2 = only overwrite originally missing calls\n"
"                     3 = only overwrite calls which are nonmissing in new file\n"
"                     4/5 = never overwrite and always overwrite, respectively\n"
"                     6 = report all mismatching calls without merging\n"
"                     7 = report mismatching nonmissing calls without merging\n"
	       );
    help_print("indiv-sort\tmerge\tbmerge\tmerge-list", &help_ctrl, 0,
"  --indiv-sort [m] : Specify family/individual ID sort order.  The following\n"
"                     three modes are currently supported:\n"
"                     * 'none'/'0' keeps individuals in the order they were\n"
"                       loaded.  This is the default for non-merge operations.\n"
"                     * 'natural'/'n' invokes \"natural sort\", e.g. 'id2' <\n"
"                       'ID3' < 'ID10'.  This is the default when merging.\n"
"                     * 'ascii'/'a' sorts in ASCII order, e.g. 'ID3' < 'id10' <\n"
"                       'id2'.\n"
"                     For now, only --make-bed and --merge/--bmerge/--merge-list\n"
"                     respect this flag.\n"
	       );
    help_print("pheno\tall-pheno", &help_ctrl, 0,
"  --pheno [fname]  : Specify alternate phenotype.\n"
"  --all-pheno      : For basic association tests, loop through all phenotypes\n"
"                     in --pheno file.\n"
	       );
    help_print("mpheno\tpheno", &help_ctrl, 0,
"  --mpheno [col]   : Specify phenotype column number in --pheno file.\n"
	       );
    help_print("pheno-name\tpheno", &help_ctrl, 0,
"  --pheno-name [c] : If phenotype file has a header row, use column with the\n"
"                     given name.\n"
	       );
    help_print("pheno-merge\tpheno", &help_ctrl, 0,
"  --pheno-merge    : If a phenotype is present in the original but not the\n"
"                     alternate file, use the original value instead of setting\n"
"                     the phenotype to missing.\n"
	       );
    help_print("covar\tcovar-name\tcovar-number", &help_ctrl, 0,
"  --covar [filename]    : Specify covariate file.\n"
"  --covar-name [names]  : Specifies covariate(s) in --covar file by name.\n"
"                          Separate multiple names with commas (spaces are not\n"
"                          permitted).  Use dashes to designate ranges.\n"
"  --covar-number [nums] : Specifies covariate(s) in --covar file by number.\n"
	       );
    help_print("within", &help_ctrl, 0,
"  --within [fname] : Specify clusters.\n"
	       );
    help_print("set\tsubset", &help_ctrl, 0,
"  --set [filename] : Specify sets.\n"
"  --subset [fname] : Specify list of subsets to extract from --set file.\n"
	       );
    help_print("prune", &help_ctrl, 0,
"  --prune          : Remove individuals with missing phenotypes.\n"
	       );
    help_print("1", &help_ctrl, 0,
"  --1              : Affection phenotypes are interpreted as 0 = unaffected,\n"
"                     1 = affected (instead of 0 = missing, 1 = unaffected,\n"
"                     2 = affected).\n"
	       );
// --map3 implicitly supported via autodetection
// --compound-genotypes automatically supported
    help_print("cow\tdog\thorse\tmouse\trice\tsheep", &help_ctrl, 0,
"  --cow/--dog/--horse/--mouse/--rice/--sheep : Specify nonhuman species.\n"
	       );
    help_print("autosome\tautosome-xy\tchr\tchr-excl", &help_ctrl, 0,
"  --autosome       : Exclude all non-autosomal markers.\n"
"  --autosome-xy    : Exclude all non-autosomal markers, except those with\n"
"                     chromosome code XY (pseudo-autosomal region of X).\n"
	       );
    help_print("chr\tchr-excl\tfrom-bp\tto-bp\tfrom-kb\tto-kb\tfrom-mb\tto-mb", &help_ctrl, 0,
"  --chr [chrs...]  : Exclude all markers not on the given chromosome(s).  Valid\n"
"                     choices for humans are 0 (unplaced), 1-22, X, Y, XY, and\n"
"                     MT.  Separate multiple chromosomes with spaces and/or\n"
"                     commas, and use a dash (no adjacent spaces permitted) to\n"
"                     denote a range, e.g. '--chr 1-4, 22, xy'.\n"
"  --chr-excl [...] : Reverse of --chr (excludes markers on listed chromosomes).\n"
	       );
    help_print("d\tsnps", &help_ctrl, 0,
"  --d [char]       : Change marker range delimiter (usually '-').\n"
	       );
    help_print("from\tto\tsnp\twindow\tfrom-bp\tto-bp\tfrom-kb\tto-kb\tfrom-mb\tto-mb", &help_ctrl, 0,
"  --from [mkr ID]  : Use ID(s) to specify a marker range to load.  When used\n"
"  --to   [mkr ID]    together, both markers must be on the same chromosome.\n"
"  --snp  [mkr ID]  : Specify a single marker to load.\n"
"  --window  [kbs]  : With --snp, loads all markers within half the specified kb\n"
"                     distance of the named marker.\n"
"  --from-bp [pos]  : Use physical position(s) to define a marker range to load.\n"
"  --to-bp   [pos]    --from-kb/--to-kb/--from-mb/--to-mb allow decimal values.\n"
"    ...              You're required to specify a single chromosome with these.\n"
	       );
    help_print("snps", &help_ctrl, 0,
"  --snps [IDs...]  : Use IDs to specify multiple marker ranges to load.  E.g.\n"
"                     '--snps rs1111-rs2222, rs3333, rs4444'.\n"
	       );
    help_print("maf", &help_ctrl, 0,
"  --maf {val}      : Minor allele frequency minimum threshold (default 0.01).\n"
"                     Note that the default threshold is only applied if --maf\n"
"                     is used without an accompanying value; if you do not\n"
"                     invoke --maf, no MAF inclusion threshold is applied.\n"
"                     Other inclusion thresholds work the same way.\n"
	       );
    help_print("max-maf\tmaf", &help_ctrl, 0,
"  --max-maf [val]  : Minor allele frequency maximum threshold.\n"
	       );
    help_print("geno\tmind", &help_ctrl, 0,
"  --geno {val}     : Maximum per-marker missing (default 0.1).\n"
"  --mind {val}     : Maximum per-person missing (default 0.1).\n"
	       );
    help_print("hwe", &help_ctrl, 0,
"  --hwe {val}      : Minimum Hardy-Weinberg disequilibrium p-value (exact),\n"
"                     default 0.001.\n"
	       );
    help_print("hwe-all\thwe", &help_ctrl, 0,
"  --hwe-all        : Given case-control data, include noncontrols in HWE test.\n"
	       );
    help_print("allow-no-sex\tmust-have-sex", &help_ctrl, 0,
"  --allow-no-sex   : Do not treat ambiguous-sex individuals as having missing\n"
"                     phenotypes in analysis commands.  (Automatic /w --no-sex.)\n"
"  --must-have-sex  : Force ambiguous-sex phenotypes to missing on --make-bed\n"
"                     and --recode.\n"
	       );
    help_print("set-hh-missing", &help_ctrl, 0,
"  --set-hh-missing : Cause --make-bed and --recode to set heterozygous haploid\n"
"                     genotypes to missing.\n"
	       );
    help_print("reference\tallele-count\tlfile", &help_ctrl, 0,
"  --reference [fn] : Specify default allele file for .lgen input.\n"
"  --allele-count   : When used with --reference, specifies that the .lgen file\n"
"                     contains reference allele counts.\n"
	       );
    help_print("with-phenotype\twrite-covar", &help_ctrl, 0,
"  --with-phenotype : Include phenotype when writing updated covariate file.\n"
	       );
    help_print("dummy-coding\twrite-covar", &help_ctrl, 0,
"  --dummy-coding   : Split categorical variables (n categories, for 2 < n < 50)\n"
"                     into n-1 binary dummy variables when writing updated\n"
"                     covariate file.\n"
	       );
    help_print("nonfounders", &help_ctrl, 0,
"  --nonfounders    : Include nonfounders in allele frequency/HWE calculations.\n"
	       );
    help_print("ppc-gap\tgenome\tZ-genome", &help_ctrl, 0,
"  --ppc-gap [val]  : Minimum number of base pairs, in thousands, between\n"
"                     informative pairs of markers used in --genome PPC test.\n"
"                     500 if unspecified.\n"
	       );
    help_print("seed", &help_ctrl, 0,
"  --seed [val...]  : Set random number seed(s).  Each value must be an integer\n"
"                     between 0 and 4294967295 inclusive.\n"
	       );
    help_print("memory", &help_ctrl, 0,
"  --memory [val]   : Size, in MB, of initial malloc attempt.  (Some operating\n"
"                     systems allow this number to exceed total physical RAM.)\n"
	       );
    help_print("threads\tthread-num", &help_ctrl, 0,
"  --threads [val]  : Maximum number of concurrent threads.\n"
	       );
    help_print("debug", &help_ctrl, 0,
"  --debug          : Enable debug logging.\n"
	       );
    help_print("extract\texclude", &help_ctrl, 0,
"  --extract [file] : Exclude all markers not in the given list.\n"
"  --exclude [file] : Exclude all markers in the given list.\n"
	       );
    help_print("keep\tremove", &help_ctrl, 0,
"  --keep [fname]   : Exclude all individuals not in the given list.\n"
"  --remove [fname] : Exclude all individuals in the given list.\n"
	       );
    help_print("maf-succ", &help_ctrl, 0,
"  --maf-succ       : Rule of succession MAF estimation (used in EIGENSTRAT).\n"
"                     Given j observations of one allele and k >= j observations\n"
"                     of the other, infer a MAF of (j+1) / (j+k+2), rather than\n"
"                     the usual j / (j+k).\n"
	       );
    help_print("ci", &help_ctrl, 0,
"  --ci [size]      : Set size of odds ratio confidence interval to report.\n"
	       );
    help_print("pfilter", &help_ctrl, 0,
"  --pfilter [val]  : Filter out association test results with higher p-values.\n"
	       );
    help_print("loop-assoc", &help_ctrl, 0,
"  --loop-assoc [f] : Run specified case/control association commands once for\n"
"                     each cluster in the file, using cluster membership as the\n"
"                     phenotype.\n"
	       );
    help_print("cell\tmodel", &help_ctrl, 0,
"  --cell [thresh]  : Specify contingency table threshold for performing all\n"
"                     --model tests.\n"
	       );
    help_print("adjust\tlambda", &help_ctrl, 0,
"  --lambda [val]   : Set genomic control lambda for --adjust.\n"
	       );
    help_print("exponent\tdistance", &help_ctrl, 0,
"  --exponent [val] : When computing genomic distances, each marker has a weight\n"
"                     of (2q(1-q))^{-val}, where q is the inferred MAF.  (Use\n"
"                     --read-freq if you want to explicitly specify some or all\n"
"                     of the MAFs.)\n"
	       );
    help_print("mperm-save\tmperm-save-all", &help_ctrl, 0,
"  --mperm-save     : Save best max(T) permutation test statistics.\n"
"  --mperm-save-all : Save all max(T) permutation test statistics.\n"
	       );
    help_print("flip\tflip-subset", &help_ctrl, 0,
"  --flip [filename]   : Flip alleles (A<->T, C<->G) for all SNP ID in the file.\n"
"  --flip-subset [fn]  : Only apply --flip to indivs in the --flip-subset file.\n"
	       );
    help_print("recode\trecode-allele", &help_ctrl, 0,
"  --recode-allele [f] : With --recode A or --recode AD, count alleles named in\n"
"                        the file (instead of the minor allele).\n"
	       );
    help_print("keep-allele-order\tmake-bed\tmerge\tbmerge\tmerge-list", &help_ctrl, 0,
"  --keep-allele-order : Keep the original allele order when creating a new\n"
"                        fileset, instead of forcing A2 to be the major allele.\n"
	       );
    help_print("merge\tbmerge\tmerge-list\tmerge-mode\tmerge-equal-pos", &help_ctrl, 0,
"  --merge-equal-pos   : Merge markers with different names but identical\n"
"                        positions.\n"
	       );
    help_print("update-chr\tupdate-cm\tupdate-map\tupdate-name", &help_ctrl, 0,
"  --update-chr  [f] {chrcol} {IDcol}  {skip} : Update marker chromosome codes.\n"
"  --update-cm   [f] {cmcol}  {IDcol}  {skip} : Update centimorgan positions.\n"
"  --update-map  [f] {bpcol}  {IDcol}  {skip} : Update marker bp positions.\n"
"  --update-name [f] {newcol} {oldcol} {skip} : Update marker IDs.\n"
	       );
    help_print("update-alleles", &help_ctrl, 0,
"  --update-alleles [fname]  : Update marker allele codes.\n"
	       );
    help_print("update-ids\tupdate-parents\tupdate-sex", &help_ctrl, 0,
"  --update-ids [fname]      : Update individual IDs, parental IDs, or sexes\n"
"  --update-parents [fname]    using the information in the provided file.  See\n"
"  --update-sex [fname]        the PLINK documentation for file format details.\n"
	       );
    help_print("simulate\tsimulate-ncases\tsimulate-ncontrols\tsimulate-prevalence", &help_ctrl, 0,
"  --simulate-ncases [num]   : Set --simulate case count (default 1000).\n"
"  --simulate-ncontrols [n]  : Set --simulate control count (default 1000).\n"
"  --simulate-prevalence [p] : Set --simulate disease prevalence (default 0.01).\n"
	       );
    help_print("simulate-qt\tsimulate-n", &help_ctrl, 0,
"  --simulate-n [num]        : Set --simulate-qt indiv count (default 1000).\n"
	       );
    help_print("simulate\tsimulate-qt\tsimulate-label\tsimulate-missing", &help_ctrl, 0,
"  --simulate-label [prefix] : Set --simulate(-qt) individual name prefix.\n"
"  --simulate-missing [freq] : Set --simulate(-qt) missing genotype frequency.\n"
	       );
    help_print("allele1234\talleleACGT\talleleacgt", &help_ctrl, 0,
"  --allele1234 <multichar>  : Interpret/recode A/C/G/T alleles as 1/2/3/4.\n"
"                              With 'multichar', converts all A/C/G/Ts in\n"
"                              allele names to 1/2/3/4s.\n"
"  --alleleACGT <multichar>  : Reverse of --allele1234.\n"
	       );
    help_print("a1-allele\treference-allele\tupdate-ref-allele\ta2-allele", &help_ctrl, 0,
"  --a1-allele [f] {a1col} {IDcol} {skip} : Force alleles in the file to A1.\n"
"  --a2-allele [f] {a2col} {IDcol} {skip} : Force alleles in the file to A2.\n"
	       );
    help_print("read-freq\tupdate-freq", &help_ctrl, 0,
"  --read-freq [filename]    : Loads MAFs from the given PLINK-style or --freqx\n"
"  --update-freq [filename]    frequency file, instead of just setting them to\n"
"                              frequencies observed in the .ped/.bed file.\n"
	       );
    help_print("parallel", &help_ctrl, 0,
"  --parallel [k] [n]        : Divide the output matrix into n pieces, and only\n"
"                              compute the kth piece.  The primary output file\n"
"                              will have the piece number included in its name,\n"
"                              e.g. wdist.rel.13 or wdist.rel.13.gz if k is 13.\n"
"                              Concatenating these files in order will yield the\n"
"                              full matrix of interest.  (Yes, this can be done\n"
"                              before unzipping.)\n"
"                              N.B. This cannot be used to directly write a\n"
"                              symmetric square matrix.  Choose the square0 or\n"
"                              triangle format instead, and postprocess as\n"
"                              necessary.\n"
	       );
    help_print("filter", &help_ctrl, 0,
"  --filter [filename] [val] : Filter individuals based on a phenotype value.\n"
	       );
    help_print("mfilter\tfilter", &help_ctrl, 0,
"  --mfilter [col]           : Specify alternate column in --filter file.\n"
	       );
    help_print("filter-cases\tfilter-controls", &help_ctrl, 0,
"  --filter-cases            : Include only cases.\n"
"  --filter-controls         : Include only controls.\n"
	       );
    help_print("filter-males\tfilter-females", &help_ctrl, 0,
"  --filter-males            : Include only males.\n"
"  --filter-females          : Include only females.\n"
	       );
    help_print("filter-founders\tfilter-nonfounders", &help_ctrl, 0,
"  --filter-founders         : Include only founders.\n"
"  --filter-nonfounders      : Include only nonfounders.\n"
	       );
    help_print("missing-genotype\tmissing-phenotype", &help_ctrl, 0,
"  --missing-genotype [char] : Code for missing genotype (normally '0').\n"
"  --missing-phenotype [val] : Numeric code for missing phenotype (normally -9).\n"
	       );
    help_print("perm-batch-size", &help_ctrl, 0,
"  --perm-batch-size [val]   : Set number of permutations per batch in QT\n"
"                              permutation tests.\n"
	       );
    help_print("output-missing-genotype\toutput-missing-phenotype", &help_ctrl, 0,
"  --output-missing-genotype [ch] : Code for missing genotype when creating new\n"
"                                   text fileset (--recode).\n"
"  --output-missing-phenotype [n] : Numeric code for missing phenotype when\n"
"                                   creating new fileset (--make-bed/--recode).\n"
	       );
    help_print("missing-code\tmissing_code\tmissing-phenotype", &help_ctrl, 0,
"  --missing-code {vals}     : Comma-separated list of missing phenotype values,\n"
"  --missing_code {vals}       for Oxford-formatted filesets (normally 'NA').\n"
	       );
    help_print("make-pheno", &help_ctrl, 0,
"  --make-pheno [file] [val] : Specify dichotomous phenotype, where cases have\n"
"                              the given value.  If the value is '*', all\n"
"                              individuals present in the phenotype file are\n"
"                              affected (and other individuals in the .ped/.fam\n"
"                              are unaffected).\n"
	       );
    help_print("tail-pheno\tgroupdist", &help_ctrl, 0,
"  --tail-pheno [Ltop] {Hbt} : Form 'low' (<= Ltop, unaffected) and 'high'\n"
"                              (greater than Hbt, affected) groups from scalar\n"
"                              phenotype data.  If Hbt is unspecified, it is set\n"
"                              equal to Ltop.  Central phenotype values are\n"
"                              treated as missing.\n"
	       );
    help_print("adjust\tgc\tlog10\tqq-plot", &help_ctrl, 0,
"  --adjust <gc> <log10> <qq-plot> : Report some multiple-testing corrections.\n"
	       );
    help_print("aperm", &help_ctrl, 0,
"  --aperm [min perms] [max perms] [alpha] [beta] [init interval] [slope] :\n"
"    This sets six parameters controlling adaptive permutation tests.  Defaults\n"
"    are 5, 1000000, 0, 0.0001, 1, and 0.001, respectively.\n"
	       );
    if (!param_ct) {
      fputs(
"\nThese flags only apply to .cnv fileset analysis:\n"
, stdout);
    }
    help_print("cnv-exclude-off-by-1", &help_ctrl, 0,
"  --cnv-exclude-off-by-1   : Exclude .cnv segments where the terminal .cnv.map\n"
"                             entry is off by 1.\n"
	       );
    help_print("cnv-del\tcnv-dup", &help_ctrl, 0,
"  --cnv-del                : Exclude all variants with more than one copy.\n"
"  --cnv-dup                : Exclude all variants with fewer than three copies.\n"
	       );
    help_print("cnv-kb\tcnv-max-kb", &help_ctrl, 0,
"  --cnv-kb [kb len]        : Exclude segments shorter than the given length.\n"
"  --cnv-max-kb [kb len]    : Exclude segments longer than the given length.\n"
	       );
    help_print("cnv-score\tcnv-max-score", &help_ctrl, 0,
"  --cnv-score [val]        : Exclude all variants with confidence score < val.\n"
"  --cnv-max-score [val]    : Exclude all variants with confidence score > val.\n"
	       );
    help_print("cnv-sites\tcnv-max-sites", &help_ctrl, 0,
"  --cnv-sites [ct]         : Exclude all segments with fewer than ct probes.\n"
"  --cnv-max-sites [ct]     : Exclude all segments with more than ct probes.\n"
	       );
    help_print("cnv-intersect\tcnv-exclude\tcnv-subset\tcnv-overlap\tcnv-region-overlap\tcnv-union-overlap\tcnv-disrupt", &help_ctrl, 0,
"  --cnv-intersect [fname]  : Only include segments which intersect a region in\n"
"                             the given region list.\n"
"  --cnv-exclude [fname]    : Exclude all segments which intersect a region in\n"
"                             the given region list.\n"
"  --cnv-count [fname]      : Specify region list for --cnv-indiv-perm\n"
"                             (optional) or --cnv-enrichment-test (required).\n"
"  --cnv-subset [fname]     : Ignore all regions in the --cnv-intersect/-exclude\n"
"                             /-count list that aren't named in the given file.\n"
"  --cnv-overlap [x]        : Only count intersections of length at least xn,\n"
"                             where n is the segment size.\n"
"  --cnv-region-overlap [x] : x >= [overlap] / [region size].\n"
"  --cnv-union-overlap [x]  : x >= [overlap] / [union size].\n"
"  --cnv-disrupt            : Only include/exclude segments with an endpoint in\n"
"                             a region.\n"
	       );
    help_print("cnv-freq-exclude-above\tcnv-freq-exclude-below\tcnv-freq-exclude-exact\tcnv-freq-include-exact\tcnv-freq-overlap\tcnv-freq-method2", &help_ctrl, 0,
"  --cnv-freq-exclude-above [k] : Exclude all segments where any portion is\n"
"                                 included by more than k total segments.\n"
"  --cnv-freq-exclude-below [k] : Exclude all segments where no portion is\n"
"                                 included by k or more total segments.\n"
"  --cnv-freq-exclude-exact [k] : Exclude all segments which have a portion\n"
"                                 included by at least k total segments, but no\n"
"                                 portion included by more.\n"
"  --cnv-freq-include-exact [k] : Reverse of --cnv-freq-exclude-exact.\n"
"  --cnv-freq-overlap {x}   : Only count portions of length at least xn, where n\n"
"                             is the segment size.\n"
"  --cnv-freq-method2 {x}   : Causes k to instead be compared against the number\n"
"                             of segments for which x >= [overlap] / [union].\n"
	       );
    help_print("cnv-test-window\tcnv-test", &help_ctrl, 0,
"  --cnv-test-window [size] : Specify window size (in kb) for CNV assoc. test.\n"
	       );
    if (!param_ct) {
      fputs(
"\nFor further documentation and support, consult the main webpage\n"
"(https://www.cog-genomics.org/wdist ) and/or the wdist-users mailing list\n"
"(https://groups.google.com/d/forum/wdist-users ).\n"
, stdout);
    }
  } while (help_ctrl.iters_left--);
  if (help_ctrl.unmatched_ct) {
    net_unmatched_ct = help_ctrl.unmatched_ct;
    printf("\nNo help entr%s for", (help_ctrl.unmatched_ct == 1)? "y" : "ies");
    col_num = (help_ctrl.unmatched_ct == 1)? 17 : 19;
    arg_uidx = 0;
    while (help_ctrl.unmatched_ct) {
      arg_uidx = next_non_set_unsafe(help_ctrl.all_match_arr, arg_uidx);
      help_ctrl.unmatched_ct--;
      if (help_ctrl.unmatched_ct) {
	if (net_unmatched_ct == 2) {
	  if (help_ctrl.param_lens[arg_uidx] + col_num > 76) {
	    printf("\n'%s'", argv[arg_uidx]);
	    col_num = 2 + help_ctrl.param_lens[arg_uidx];
	  } else {
	    printf(" '%s'", argv[arg_uidx]);
	    col_num += 3 + help_ctrl.param_lens[arg_uidx];
	  }
	} else {
	  if (help_ctrl.param_lens[arg_uidx] + col_num > 75) {
	    printf("\n'%s',", argv[arg_uidx]);
	    col_num = 3 + help_ctrl.param_lens[arg_uidx];
	  } else {
	    printf(" '%s',", argv[arg_uidx]);
	    col_num += 4 + help_ctrl.param_lens[arg_uidx];
	  }
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
	if (help_ctrl.param_lens[arg_uidx] + col_num > 75) {
	  printf("\n'%s'.\n", argv[arg_uidx]);
	} else {
	  printf(" '%s'.\n", argv[arg_uidx]);
	}
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

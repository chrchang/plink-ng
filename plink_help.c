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

#define MAX_EQUAL_HELP_PARAMS 22

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
	    arg_uidx = next_unset_unsafe(help_ctrl_ptr->all_match_arr, arg_uidx);
	    for (cur_param_idx = 0; cur_param_idx < cur_param_ct; cur_param_idx++) {
	      if (!strcmp(cur_param_start[cur_param_idx], help_ctrl_ptr->argv[arg_uidx])) {
		SET_BIT(help_ctrl_ptr->perfect_match_arr, arg_uidx);
		SET_BIT(help_ctrl_ptr->prefix_match_arr, arg_uidx);
		SET_BIT(help_ctrl_ptr->all_match_arr, arg_uidx);
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
		  SET_BIT(help_ctrl_ptr->prefix_match_arr, arg_uidx);
		  SET_BIT(help_ctrl_ptr->all_match_arr, arg_uidx);
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
		SET_BIT(help_ctrl_ptr->all_match_arr, arg_uidx);
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
"Most " PROG_NAME_CAPS " runs require exactly one main input fileset.  The following flags\n"
"are available for defining its form and location:\n\n"
, stdout);
  }
  do {
    help_print("bfile\tbed\tbim\tfam", &help_ctrl, 1,
"  --bfile {prefix} : Specify .bed + .bim + .fam prefix (default '" PROG_NAME_STR "').\n"
"  --bed [filename] : Specify full name of .bed file.\n"
"  --bim [filename] : Specify full name of .bim file.\n"
"  --fam [filename] : Specify full name of .fam file.\n\n"
	       );
    help_print("file\tped\tmap", &help_ctrl, 1,
"  --file {prefix}  : Specify .ped + .map filename prefix (default '" PROG_NAME_STR "').\n"
"  --ped [filename] : Specify full name of .ped file.\n"
"  --map [filename] : Specify full name of .map file.\n\n"
	       );
    help_print("bfile\tfam\tfile\tped\tno-fid\tno-parents\tno-sex\tno-pheno", &help_ctrl, 1,
"  --no-fid         : .fam/.ped file does not contain column 1 (family ID).\n"
"  --no-parents     : .fam/.ped file does not contain columns 3-4 (parents).\n"
"  --no-sex         : .fam/.ped file does not contain column 5 (sex).\n"
"  --no-pheno       : .fam/.ped file does not contain column 6 (phenotype).\n\n"
	       );
    help_print("tfile\ttped\ttfam", &help_ctrl, 1,
"  --tfile {prefix} : Specify .tped + .tfam filename prefix (default '" PROG_NAME_STR "').\n"
"  --tped [fname]   : Specify full name of .tped file.\n"
"  --tfam [fname]   : Specify full name of .tfam file.\n\n"
	       );
    help_print("lfile\treference\tallele-count", &help_ctrl, 1,
"  --lfile {prefix} : Specify .lgen + .map + .fam (long-format fileset) prefix.\n"
"  --reference [fn] : Specify default allele file accompanying --lfile input.\n"
"  --allele-count   : When used with --lfile + --reference, specifies that the\n"
"                     .lgen file contains reference allele counts.\n\n"
	       );
    help_print("vcf\tbcf", &help_ctrl, 1,
"  --vcf [filename] : Specify full name of .vcf or .vcf.gz file.\n"
"  --bcf [filename] : Specify full name of .bcf file.\n\n"
	       );
    help_print("file\ttfile\tlfile\t23file\tkeep-autoconv", &help_ctrl, 1,
"  --keep-autoconv  : With --file/--tfile/--lfile/--vcf/--23file, don't delete\n"
"                     autogenerated binary fileset at end of run.\n\n"
	       );
    help_print("23file", &help_ctrl, 1,
"  --23file [fname] {family ID} {indiv. ID} {sex} {pheno} {pat. ID} {mat. ID} :\n"
"    Specify 23andMe input file.\n\n"
	       );
    /*
    help_print("data\tgen\tsample", &help_ctrl, 1,
"  --data {prefix}  : Specify Oxford .gen + .sample prefix (default '" PROG_NAME_STR "').\n"
"  --gen [filename] : Specify full name of .gen file.\n"
"  --sample [fname] : Specify full name of .sample file.\n\n"
    	       );
    */
#ifndef STABLE_BUILD
    help_print("cfile\tcnv-list\tgfile", &help_ctrl, 1,
"  --cfile [prefix] : Specify .cnv + .fam + .cnv.map (segmental CNV) prefix.\n"
"  --cnv-list [fn]  : Specify full name of .cnv file.\n"
"  --gfile [prefix] : Specify .gvar + .fam + .map (genetic variant) prefix.\n\n"
	       );
#endif
    help_print("grm\tgrm-bin\trel-cutoff\tgrm-cutoff", &help_ctrl, 1,
"  --grm {prefix}   : Specify .grm.gz + .grm.id (GCTA rel. matrix) prefix.\n"
"  --grm-bin {prfx} : Specify .grm.bin + .grm.N.bin + .grm.id (GCTA triangular\n"
"                     binary relationship matrix) filename prefix.\n\n"
	       );
    help_print("dummy", &help_ctrl, 1,
"  --dummy [indiv ct] [SNP ct] {missing geno freq} {missing pheno freq}\n"
"          <acgt | 1234 | 12> <scalar-pheno>\n"
"    This generates a fake input dataset with the specified number of\n"
"    individuals and SNPs.  By default, the missing genotype and phenotype\n"
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
"Output files have names of the form '" PROG_NAME_STR ".{extension}' by default.  You can\n"
"change the '" PROG_NAME_STR "' prefix with\n\n"
, stdout);
    }
    help_print("out", &help_ctrl, 1,
"  --out [prefix]   : Specify prefix for output files.\n\n"
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
    help_print("recode\trecode12\ttab\ttranspose\trecode-lgen\trecodeAD\trecodead\trecodeA\trecodea\trecode-rlist\trecode-allele\tlist\twith-reference\trecode-vcf\tfid\tiid\trecode-beagle\trecode-bimbam\trecode-fastphase\trecodeHV\trecodehv\trecode-structure", &help_ctrl, 1,
"  --recode <12> <compound-genotypes> <23 | A | AD | beagle | bimbam |\n"
"           bimbam-1chr | fastphase | fastphase-1chr | HV | HV-1chr | lgen |\n"
"           lgen-ref | list | oxford | rlist | structure | transpose | vcf |\n"
"           vcf-fid | vcf-iid> <tab | tabx | spacex>\n"
"    Create a new text fileset with all filters applied.  By default, the\n"
"    fileset consists of a .ped and a .map file, readable with --file.\n"
"    * The '12' modifier causes all alleles to be coded as 1s and 2s.\n"
"    * The 'compound-genotypes' modifier removes the space between pairs of\n"
"      genotype codes for the same variant.\n"
"    * The '23' modifier causes a 23andMe-formatted file to be generated.  This\n"
"      can only be used on a single individual's data (--keep may be handy).\n"
"    * The 'AD' modifier causes an additive + dominant component file, suitable\n"
"      for loading from R, to be generated.  If you don't want the dominant\n"
"      component, use 'A' instead.\n"
"    * The 'beagle' modifier causes per-autosome .dat and .map files, readable\n"
"      by BEAGLE, to be generated.\n"
"    * The 'bimbam' modifier causes a BIMBAM-formatted fileset to be generated.\n"
"      If your input data only contains one chromosome, you can use\n"
"      'bimbam-1chr' instead to write a two-column .pos.txt file.\n"
"    * The 'fastphase' modifier causes per-chromosome fastPHASE files to be\n"
"      generated.  If your input data only contains one chromosome, you can use\n"
"      'fastphase-1chr' instead to exclude the chromosome number from the file\n"
"      extension.\n"
"    * The 'HV' modifier causes a Haploview-format .ped + .info fileset to be\n"
"      generated per chromosome.  'HV-1chr' is analogous to 'fastphase-1chr'.\n"
"    * The 'lgen' modifier causes a long-format fileset (loadable with --lfile)\n"
"      to be generated, while 'lgen-ref' generates a (usually) smaller\n"
"      long-format fileset loadable with --lfile + --reference.\n"
"    * The 'list' modifier creates a genotype-based list, while 'rlist'\n"
"      creates a rare-genotype fileset.\n"
"    * 'oxford' causes an Oxford-format .gen + .sample fileset to be generated.\n"
"    * The 'structure' modifier causes a Structure-format file to be generated.\n"
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
"    If a --covar file is loaded, --make-bed and --recode automatically generate\n"
"    an updated version (with all filters applied).  However, if you do not wish\n"
"    to simultaneously generate a new genotype file, you can use --write-covar\n"
"    to just produce a pruned covariate file.\n\n"
	       );
    help_print("write-cluster", &help_ctrl, 1,
"  --write-cluster <omit-unassigned>\n"
"    If a --within file is loaded, this generates another cluster file (with all\n"
"    filters applied).  The 'omit-unassigned' modifier causes unclustered\n"
"    individuals to be omitted from the file; otherwise their cluster is \"NA\".\n\n"
	       );
#ifndef STABLE_BUILD
    help_print("write-set\tset-table", &help_ctrl, 1,
"  --write-set\n"
"  --set-table\n"
"    If sets have been defined, --write-set dumps 'END'-terminated set\n"
"    membership lists to {output prefix}.set, while --set-table writes a\n"
"    variant-by-set membership table to {output prefix}.set.table.\n\n"
	       );
#endif
    help_print("merge\tbmerge\tmerge-list\tmerge-mode", &help_ctrl, 1,
"  --merge [.ped filename] [.map filename]\n"
"  --merge [text fileset prefix]\n"
"  --bmerge [.bed filename] [.bim filename] [.fam filename]\n"
"  --bmerge [binary fileset prefix]\n"
"    Merge the given fileset with the initially loaded fileset, writing the\n"
"    result to {output prefix}.bed + .bim + .fam.  (It is no longer necessary to\n"
"    simultaneously specify --make-bed.\n"
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
    help_print("write-snplist\tlist-23-indels", &help_ctrl, 1,
"  --write-snplist\n"
"  --list-23-indels\n"
"    --write-snplist writes a .snplist file listing the names of all variants\n"
"    which pass the filters and inclusion thresholds you've specified, while\n"
"    --list-23-indels writes the subset with 23andMe-style indel calls (I/D\n"
"    allele codes).\n\n"
	       );
    help_print("freq\tfreqx\tfrqx\tcounts", &help_ctrl, 1,
"  --freq <counts>\n"
"  --freqx\n"
"    --freq generates a basic allele frequency (or count, if the 'counts'\n"
"    modifier is present) report.  This can be combined with --within to produce\n"
"    a cluster-stratified allele frequency/count report instead.\n"
"    --freqx generates a more detailed genotype count report, designed for use\n"
"    with --read-freq.\n\n"
		);
    help_print("missing", &help_ctrl, 1,
"  --missing\n"
"    Generate individual- and variant-based missing data reports.  If clusters\n"
"    are defined, the variant-based report is cluster-stratified.\n\n"
	       );
    help_print("hardy", &help_ctrl, 1,
"  --hardy\n"
"    Generate a Hardy-Weinberg exact test p-value report.  (This does NOT\n"
"    simultaneously filter on the p-value any more; use --hwe for that.)\n\n"
	       );
    help_print("ibc\thet", &help_ctrl, 1,
"  --ibc\n"
"    Calculate inbreeding coefficients in three different ways.  (The second is\n"
"    the excess homozygosity coefficient calculated by the old --het command.)\n"
"    * For more details, see Yang J, Lee SH, Goddard ME and Visscher PM.  GCTA:\n"
"      a tool for Genome-wide Complex Trait Analysis.  Am J Hum Genet. 2011 Jan\n"
"      88(1): 76-82.  This paper also describes the relationship matrix\n"
"      computation we implement.\n\n"
	       );
    help_print("indep\tindep-pairwise", &help_ctrl, 1,
"  --indep [window size]<kb> [step size (site ct)] [VIF threshold]\n"
"  --indep-pairwise [window size]<kb> [step size (site ct)] [r^2 threshold]\n"
"    Generate a list of markers in approximate linkage equilibrium.  With the\n"
"    'kb' modifier, the window size is in kilobase instead of site count units.\n"
"    (Pre-'kb' space is optional, i.e. '--indep-pairwise 500 kb 5 0.5' and\n"
"    '--indep-pairwise 500kb 5 0.5' have the same effect.)\n"
"    Note that you need to rerun " PROG_NAME_CAPS " using --extract or --exclude on the\n"
"    .prune.in/.prune.out file to apply the list to another computation.\n\n"
		);
    help_print("r\tr2\tmatrix", &help_ctrl, 1,
"  --r <square | square0 | triangle | inter-chr> <gz | bin> <single-prec>\n"
"      <spaces> <yes-really>\n"
"  --r2 <square | square0 | triangle | inter-chr> <gz | bin> <single-prec>\n"
"       <spaces> <yes-really>\n"
"    LD statistic reports.  --r yields raw inter-variant correlations, while\n"
"    --r2 reports their squares.  You can request results for all pairs in\n"
"    matrix format (if you specify 'bin' or one of the shape modifiers), all\n"
"    pairs in table format ('inter-chr'), or a limited window in table format\n"
"    (default).\n"
"    * The 'gz' modifier causes the output text file to be gzipped.\n"
"    * 'bin' causes the output matrix to be written in binary format.  The\n"
"      matrix is square if no shape is explicitly specified.\n"
"    * 'single-prec' causes the computation to use (and with 'bin', results to\n"
"      be saved as) single-precision instead of double-precision floating point\n"
"      numbers.\n"
"    * By default, text matrices are tab-delimited; 'spaces' switches this.\n"
"    * Since the resulting file can easily be huge, you're required to add the\n"
"      'yes-really' modifier when requesting an unfiltered, non-distributed all\n"
"      pairs computation on more than 400k variants.\n"
"    * These computations can be subdivided with --parallel (even when the\n"
"      'square' modifier is active).\n\n"
	       );
    help_print("distance", &help_ctrl, 1,
"  --distance <square | square0 | triangle> <gz | bin> <ibs> <1-ibs> <allele-ct>\n"
"             <flat-missing>\n"
"    Write a lower-triangular tab-delimited table of (weighted) genomic\n"
"    distances in allele count units to {output prefix}.dist, and a list of the\n"
"    corresponding family/individual IDs to {output prefix}.dist.id.  The first\n"
"    row of the .dist file contains a single {genome 1-genome 2} distance, the\n"
"    second row has the {genome 1-genome 3} and {genome 2-genome 3} distances in\n"
"    that order, etc.\n"
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
"    These deprecated commands are equivalent to '--distance 1-ibs flat-missing\n"
"    square' and '--distance ibs flat-missing square', respectively, except that\n"
"    they generate space- instead of tab-delimited text matrices.\n\n"
		);
    help_print("make-rel", &help_ctrl, 1,
"  --make-rel <square | square0 | triangle> <gz | bin> <cov | ibc2 | ibc3>\n"
"             <single-prec>\n"
"    Write a lower-triangular variance-standardized relationship (coancestry)\n"
"    matrix to {output prefix}.rel, and corresponding IDs to\n"
"    {output prefix}.rel.id.\n"
"    * 'square', 'square0', 'triangle', 'gz', and 'bin' act as they do on\n"
"      --distance.\n"
"    * The 'cov' modifier removes the variance standardization step, causing a\n"
"      covariance matrix to be calculated instead.\n"
"    * By default, the diagonal elements in the relationship matrix are based on\n"
"      --ibc's Fhat1; use the 'ibc2' or 'ibc3' modifiers to base them on Fhat2\n"
"      or Fhat3 instead.\n"
"    * " PROG_NAME_CAPS " normally performs this calculation with double-precision floating\n"
"      point numbers.  The 'single-prec' modifier switches to single-precision\n"
"      arithmetic, which sacrifices a bit of accuracy to decrease memory usage\n"
"      and reduce computation time.\n"
"    * The computation can be subdivided with --parallel.\n"
               );
    help_print("make-grm\tmake-grm-bin\tgrm\tgrm-bin\tmake-grm-gz", &help_ctrl, 1,
"  --make-grm-gz <no-gz> <cov | ibc2 | ibc3> <single-prec>\n"
"  --make-grm-bin <cov | ibc2 | ibc3>\n"
"    --make-grm-gz writes the relationships in GCTA's original gzipped list\n"
"    format, which describes one pair per line, while --make-grm-bin writes them\n"
"    in GCTA 1.1+'s single-precision triangular binary format.  Note that these\n"
"    formats explicitly report the number of valid observations (where neither\n"
"    individual has a missing call) for each pair, which is useful input for\n"
"    some scripts.\n"
"    These computations can be subdivided with --parallel.\n\n"
	       );
    help_print("rel-cutoff\tgrm-cutoff", &help_ctrl, 1,
"  --rel-cutoff {val}\n"
"    (alias: --grm-cutoff)\n"
"    Exclude one member of each pair of individuals with relatedness greater\n"
"    than the given cutoff value (default 0.025).  If no later operation will\n"
"    cause the list of remaining individuals to be written to disk, this will\n"
"    save it to {output prefix}.rel.id.\n"
"    Note that maximizing the remaining sample size is equivalent to the NP-hard\n"
"    maximum independent set problem, so we use a greedy algorithm instead of\n"
"    guaranteeing optimality.  (Use the --make-rel and --keep/--remove flags if\n"
"    you want to try to do better.)\n\n"
	       );
    help_print("ibs-test\tgroupdist", &help_ctrl, 1,
"  --ibs-test {permutation count}\n"
"  --groupdist {iters} {d}\n"
"    Given case/control phenotype data, these commands consider three subsets of\n"
"    the distance matrix: pairs of affected individuals, affected-unaffected\n"
"    pairs, and pairs of unaffected individuals.  Each of these subsets has a\n"
"    distribution of pairwise genomic distances; --ibs-test uses permutation to\n"
"    estimate p-values re: which types of pairs are most similar, while\n"
"    --groupdist focuses on the differences between the centers of these\n"
"    distributions and estimates standard errors via delete-d jackknife.\n\n"
	       );
    help_print("regress-distance\tregress-rel", &help_ctrl, 1,
"  --regress-distance {iters} {d}\n"
"    Linear regression of pairwise genomic distances on pairwise average\n"
"    phenotypes and vice versa, using delete-d jackknife for standard errors.  A\n"
"    scalar phenotype is required.\n"
"    * With less than two parameters, d is set to {number of people}^0.6 rounded\n"
"      down.  With no parameters, 100k iterations are run.\n"
"  --regress-rel {iters} {d}\n"
"    Linear regression of pairwise genomic relationships on pairwise average\n"
"    phenotypes, and vice versa.  Defaults for iters and d are the same as for\n"
"    --regress-distance.\n\n"
	       );
    help_print("genome\tZ-genome\trel-check\timpossible\tnudge\tgenome-full\tunbounded", &help_ctrl, 1,
"  --genome <gz> <rel-check> <full> <unbounded> <nudge>\n"
"    Generate an identity-by-descent report.\n"
"    * The 'rel-check' modifier excludes pairs of individuals with different\n"
"      FIDs from the final report.\n"
"    * 'full' adds raw pairwise comparison data to the report.\n"
"    * The P(IBD=0/1/2) estimator employed by this command sometimes yields\n"
"      numbers outside the range [0,1]; by default, these are clipped.  The\n"
"      'unbounded' modifier turns off this clipping.\n"
"    * Then, when PI_HAT^2 < P(IBD=2), 'nudge' adjusts the final P(IBD=0/1/2)\n"
"      estimates to a theoretically possible configuration.\n"
"    * The computation can be subdivided with --parallel.\n\n"
		);
    help_print("homozyg\thomozyg-snp\thomozyg-kb\thomozyg-density\thomozyg-gap\thomozyg-het\thomozyg-window-snp\thomozyg-window-het\thomozyg-window-missing\thomozyg-window-threshold", &help_ctrl, 1,
"  --homozyg <group | group-verbose> <consensus-match> <extend>\n"
"            <subtract-1-from-lengths>\n"
"  --homozyg-snp [min var count]\n"
"  --homozyg-kb [min length]\n"
"  --homozyg-density [max inverse density (kb/var)]\n"
"  --homozyg-gap [max internal gap kb length]\n"
"  --homozyg-het [max hets]\n"
"  --homozyg-window-snp [scanning window size]\n"
"  --homozyg-window-het [max hets in scanning window hit]\n"
"  --homozyg-window-missing [max missing calls in scanning window hit]\n"
"  --homozyg-window-threshold [min scanning window hit rate]\n"
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
"      * By default, segment bp lengths are calculated as [end bp position] -\n"
"        [start bp position] + 1.  Therefore, reports normally differ slightly\n"
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
"  --cluster <cc> <group-avg | old-tiebreaks> <missing> <only2>\n"
"    Cluster individuals using a pairwise similarity statistic (normally IBS).\n"
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
    help_print("neighbour\tneighbor", &help_ctrl, 1,
"  --neighbour [n1] [n2]\n"
"    (alias: --neighbor)\n"
"    Report IBS distances from each individual to their n1th- to n2th-nearest\n"
"    neighbors, associated Z-scores, and the identities of those neighbors.\n"
"    Useful for outlier detection.\n\n"
	       );
    help_print("assoc\tmodel\tfisher\tperm\tmperm\tperm-count\tcounts\tp2\tmodel-dom\tmodel-gen\tmodel-rec\tmodel-trend\tgenedrop\tqt-means\ttrend", &help_ctrl, 1,
"  --assoc <perm | mperm=[value]> <perm-count> <fisher> <counts>\n"
	       /*
"  --assoc <perm | mperm=[value]> <genedrop> <perm-count> <fisher> <counts> <p2>\n"
	       */
"  --assoc <perm | mperm=[value]> <perm-count> <qt-means> <lin>\n"
"  --model <perm | mperm=[value]> <perm-count> <fisher | trend-only>\n"
      /*
"  --model <perm | mperm=[value]> <genedrop> <perm-count> <fisher | trend-only>\n"
										  */
"          <dom | rec | gen | trend>\n"
"    Basic association analysis report.\n"
"    Given a case/control phenotype, --assoc performs a 1df chi-square allelic\n"
"    test, while --model performs 4 other tests as well (1df dominant gene\n"
"    action, 1df recessive gene action, 2df genotypic, Cochran-Armitage trend).\n"
"    * With 'fisher', Fisher's exact test is used to generate p-values.\n"
"    * 'perm' causes an adaptive permutation test to be performed.\n"
"    * 'mperm=[value]' causes a max(T) permutation test with the specified\n"
"      number of replications to be performed.\n"
	       /*
"    * 'genedrop' causes offspring genotypes to be regenerated via gene-dropping\n"
"      in the permutation test.\n"
	       */
"    * 'perm-count' causes the permutation test report to include counts instead\n"
"      of frequencies.\n"
"    * 'counts' causes --assoc to report allele counts instead of frequencies.\n"
	       /*
"    * 'p2' changes the --assoc permutation test used (see PLINK documentation).\n"
	       */
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
"    Given both a quantitative phenotype and a case/control covariate loaded\n"
"    with --covar defining two groups, --gxe compares the regression coefficient\n"
"    derived from considering only members of one group to the regression\n"
"    coefficient derived from considering only members of the other.  By\n"
"    default, the first covariate in the --covar file defines the groups; use\n"
"    e.g. '--gxe 3' to base them on the third covariate instead.\n\n"
	       );
    help_print("linear\tlogistic\tperm\tmperm\tperm-count\tgenotypic\thethom\tdominant\trecessive\tno-snp\thide-covar\tsex\tno-x-sex\tinteraction\tstandard-beta\tbeta", &help_ctrl, 1,
               /*
"  --linear <perm | mperm=[value]> <genedrop> <perm-count>\n"
               */
#ifndef NOLAPACK
"  --linear <perm | mperm=[value]> <perm-count>\n"
"           <genotypic | hethom | dominant | recessive | no-snp> <hide-covar>\n"
"           <sex | no-x-sex> <interaction> <beta> <standard-beta>\n"
#endif
	       /*
"  --logistic <perm | mperm=[value]> <genedrop> <perm-count>\n"
	       */
"  --logistic <perm | mperm=[value]> <perm-count>\n"
"             <genotypic | hethom | dominant | recessive | no-snp> <hide-covar>\n"
"             <sex | no-x-sex> <interaction> <beta>\n"
"    Multi-covariate association analysis on a quantitative (--linear) or\n"
"    case/control (--logistic) phenotype.  Normally used with --covar.\n"
"    * 'perm' normally causes an adaptive permutation test to be performed on\n"
"      the main effect, while 'mperm=[value]' starts a max(T) permutation test.\n"
	       /*
"    * 'genedrop' causes offspring genotypes to be regenerated via gene-dropping\n"
"      in the permutation test.\n"
	       */
"    * 'perm-count' causes the permutation test report to include counts instead\n"
"      of frequencies.\n"
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
"    * For logistic regressions, the 'beta' modifier causes regression\n"
"      coefficients instead of odds ratios to be reported.\n"
"    * With --linear, the 'standard-beta' modifier standardizes the phenotype\n"
"      and all predictors to zero mean and unit variance before regression.\n\n"
	       );
    help_print("lasso", &help_ctrl, 1,
"  --lasso [h2 estimate] {min lambda} <report-zeroes>\n"
"    Estimate variant effect sizes via the LASSO regression discussed in\n"
"    Vattikuti S, Lee J, Hsu S, Chow CC (2013) Application of compressed sensing\n"
"    to genome wide association studies and genomic selection\n"
"    (http://arxiv.org/abs/1310.2264 ).  You must provide an additive\n"
"    heritability estimate to calibrate the regression.\n"
"    Note that this method may require a very large sample size (e.g. hundreds\n"
"    of thousands) to be effective on complex polygenic traits.\n\n"
	       );
#ifndef STABLE_BUILD
#ifndef NOLAPACK
    help_print("unrelated-heritability", &help_ctrl, 1,
"  --unrelated-heritability <strict> {tol} {initial covg} {initial covr}\n"
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
"      Metabolic Syndrome Traits.  PLoS Genet 8(3): e1002637.\n"
"      doi:10.1371/journal.pgen.1002637\n\n"
	       );
#endif
#endif
    help_print("fast-epistasis\tepistasis\tset-test\tset-by-all\tcase-only\tnop\tepistasis-summary-merge", &help_ctrl, 1,
"  --fast-epistasis <no-ueki | joint-effects> <case-only>\n"
"                   <set-by-set | set-by-all> <nop>\n"
"  --epistasis <set-by-set | set-by-all>\n"
"    Scan for epistatic interactions.  --fast-epistasis inspects 3x3 joint\n"
"    genotype count tables and only applies to case/control phenotypes, while\n"
"    --epistasis performs linear or logistic regression.\n"
"    * --fast-epistasis normally applies the variance and empty cell corrections\n"
"      described in Ueki M, Cordell HJ (2012) Improved statistics for\n"
"      genome-wide interaction analysis.  To disable them, add 'no-ueki'.\n"
"    * joint-effects' entirely replaces the original fast-epistasis statistic\n"
"      with the Ueki-Cordell 'joint effects' test.\n"
"    * 'case-only' requests a case-only instead of a case/control test.\n"
"    * By default, all pairs of variants across the entire genome are tested.\n"
"      To just test pairs of variants within a single set, add the 'set-by-set'\n"
"      modifier and load exactly one set with --set/--make-set; with exactly two\n"
"      sets loaded, all variants in one set are tested against all variants in\n"
"      the other.  'set-by-all' tests all variants in one set against the entire\n"
"      genome instead.\n"
"    * 'nop' strips p-values from the main report.\n"
"    * These computations can be subdivided with --parallel; however...\n"
"  --epistasis-summary-merge [common file prefix] [ct]\n"
"    When a --[fast-]epistasis job is subdivided with --parallel, the main\n"
"    report can be assembled at the end by applying Unix 'cat' in the usual\n"
"    manner, but the .summary.1, .summary.2, ... files require a specialized\n"
"    merge.  --epistasis-summary-merge takes care of the latter.\n\n"
	       );
    help_print("twolocus", &help_ctrl, 1,
"  --twolocus [variant ID] [variant ID]\n"
"    Two-locus genotype count report.\n\n"
	       );
    /*
    help_print("regress-pcs\tregress-pcs-distance", &help_ctrl, 1,
"  --regress-pcs [.evec or .eigenvec filename] <normalize-pheno> <sex-specific>\n"
"                <clip> {max PCs}\n"
"    Linear regression of phenotypes and genotypes on the given list of\n"
"    principal components (produced by SMARTPCA or GCTA).  Output is currently a\n"
"    .gen + .sample fileset in the Oxford IMPUTE/SNPTEST v2 format.\n"
"    * The 'normalize-pheno' modifier converts all phenotype residuals to\n"
"      Z-scores.  When combined with 'sex-specific', the Z-scores are evaluated\n"
"      separately by sex.\n"
"    * The 'clip' modifier clips out-of-range genotype residuals.  Without it,\n"
"      they are represented as negative probabilities in the .gen file, which\n"
"      are invalid input for some programs.\n"
"    * By default, principal components beyond the 20th are ignored; change this\n"
"      by setting the max PCs parameter.\n\n"
      );
    help_print("regress-pcs\tdistance\tregress-pcs-distance", &help_ctrl, 1,
"  --regress-pcs-distance [.evec/.eigenvec file] <normalize-pheno>\n"
"                         <sex-specific> {max PCs} <square | square0 | triangle>\n"
"                         <gz | bin> <ibs> <1-ibs> <allele-ct> <flat-missing>\n"
"    High-speed combination of --regress-pcs and --distance (no .gen + .sample\n"
"    fileset is written to disk).\n\n"
	       );
				     */
#ifndef STABLE_BUILD
    help_print("cnv-make-map", &help_ctrl, 1,
"  --cnv-make-map <short>\n"
"    Given a .cnv file, this generates the corresponding .cnv.map file needed\n"
"    by " PROG_NAME_CAPS "'s other CNV analysis commands.  The 'short' modifier causes\n"
"    causes entries needed by old PLINK versions to be omitted.  (This is now\n"
"    automatically invoked, with 'short', when necessary.)\n\n"
	       );
    help_print("cnv-write", &help_ctrl, 1,
"  --cnv-write <freq>\n"
"    Write a new .cnv fileset, after applying all requested filters.  The 'freq'\n"
"    modifier (which must be used with --cnv-freq-method2) causes an additional\n"
"    'FREQ' field to be written with CNV-CNV overlap counts.\n\n"
	       );
    /*
    help_print("cnv-check-no-overlap", &help_ctrl, 1,
"  --cnv-check-no-overlap\n"
"    Given a .cnv fileset, this checks for within-individual CNV overlaps.\n\n"
	       );
    */
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
#endif
    if (!param_ct) {
      fputs(
"The following other flags are supported.  (Order of operations is described at\n"
#ifdef PLINK_BUILD
"[website TBD].)\n"
#else
"https://www.cog-genomics.org/wdist/order .)\n"
#endif
, stdout);
    }
    help_print("script\trerun", &help_ctrl, 0,
"  --script [fname] : Include command-line options from file.\n"
"  --rerun {log}    : Rerun commands in log (default '" PROG_NAME_STR ".log').\n"
	       );
    help_print("silent", &help_ctrl, 0,
"  --silent         : Suppress output to console.\n"
	       );
    help_print("missing-genotype", &help_ctrl, 0,
"  --missing-genotype [char] : Set missing genotype code (normally '0').\n"
	       );
    help_print("vcf\tbcf\tdouble-id\tconst-fid\tid-delim", &help_ctrl, 0,
"  --double-id      : Set both family and individual IDs to the VCF sample ID.\n"
"  --const-fid {ID} : Set all FIDs to the given constant (default '0').\n"
"  --id-delim {d}   : Parse sample IDs as [FID][d][IID] (default delim '_').\n"
	       );
    help_print("vcf\tbcf\tbiallelic-only\tvcf-min-qual\tvcf-filter", &help_ctrl, 0,
"  --biallelic-only <strict> <list> : Skip VCF variants with 2+ alt. alleles.\n"
"  --vcf-min-qual [val]             : Skip VCF variants with low/missing QUAL.\n"
"  --vcf-filter {exception(s)...}   : Skip variants which have FILTER failures.\n"
	       );
    help_print("allow-extra-chr", &help_ctrl, 0,
"  --allow-extra-chr <0>     : Permit unrecognized chromosome codes.  The '0'\n"
"                              modifier causes them to be treated as if they had\n"
"                              been set to zero.\n"
               );
    help_print("chr-set\tcow\tdog\thorse\tmouse\trice\tsheep\tautosome-num", &help_ctrl, 0,
"  --chr-set [autosome ct] <no-x> <no-y> <no-xy> <no-mt> :\n"
"    Specify a nonhuman chromosome set.  The first parameter sets the number of\n"
"    diploid autosome pairs if positive, or haploid chromosomes if negative.\n"
"    Given diploid autosomes, the remaining modifiers indicate the absence of\n"
"    the named non-autosomal chromosomes.\n"
"  --cow/--dog/--horse/--mouse/--rice/--sheep : Shortcuts for those species.\n"
"  --autosome-num [value]    : Alias for '--chr-set [value] no-y no-xy no-mt'.\n"
	       );
    help_print("23file\t23file-convert-xy\t23file-make-xylist", &help_ctrl, 0,
"  --23file-convert-xy {f}   : Separate out XY pseudo-autosomal region.  A\n"
"                              variant list (from e.g. --23file-make-xylist) is\n"
"                              necessary to use this on a female genome.\n"
"  --23file-make-xylist      : Given a male 23andMe genome, list XY pseudo-\n"
"                              autosomal region vars in {output prefix}.xylist.\n"
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
    /*
    help_print("missing-code\tmissing_code\tmissing-phenotype", &help_ctrl, 0,
"  --missing-code {vals}     : Comma-separated list of missing phenotype values\n"
"    (alias: --missing_code)   for Oxford-formatted filesets (normally 'NA').\n"
	       );
    */
    help_print("pheno\tall-pheno\tmpheno\tpheno-name\tpheno-merge", &help_ctrl, 0,
"  --pheno [fname]  : Load phenotype data from the specified file, instead of\n"
"                     using the values in the main input fileset.\n"
"  --all-pheno      : For basic association tests, loop through all phenotypes\n"
"                     in --pheno file.\n"
"  --mpheno [col]   : Specify phenotype column number in --pheno file.\n"
"  --pheno-name [c] : If --pheno file has a header row, use column with the\n"
"                     given name.\n"
"  --pheno-merge    : When the main input fileset contains an phenotype value\n"
"                     for an individual, but the --pheno file does not, use the\n"
"                     original value instead of treating the phenotype as\n"
"                     missing.\n"
	       );
    help_print("missing-phenotype\t1", &help_ctrl, 0,
"  --missing-phenotype [v] : Set missing phenotype value (normally -9).\n"
"  --1                     : Expect case/control phenotypes to be coded as\n"
"                            0 = control, 1 = case, instead of the usual\n"
"                            0 = missing, 1 = control, 2 = case.\n"
	       );
    help_print("make-pheno\tpheno", &help_ctrl, 0,
"  --make-pheno [fn] [val] : Define a new case/control phenotype.  If the val\n"
"                            parameter is '*', all individuals listed in the\n"
"                            the given file are cases, and everyone else is a\n"
"                            control.  Otherwise, only individuals both listed\n"
"                            in the file and with third column entry identical\n"
"                            to the val parameter are cases.\n"
	       );
    help_print("tail-pheno\tgroupdist\tpheno", &help_ctrl, 0,
"  --tail-pheno [Lt] {Hbt} : Downcode a scalar phenotype to a case/control\n"
"                            phenotype.  All individuals with phenotype values\n"
"                            greater than Hbt are cases, and all with values\n"
"                            less than or equal to Lt are controls.  If Hbt is\n"
"                            unspecified, it is equal to Lt; otherwise,\n"
"                            in-between phenotype values are set to missing.\n"
	       );
    help_print("covar\tcovar-name\tcovar-number", &help_ctrl, 0,
"  --covar [filename] <keep-pheno-on-missing-cov> : Specify covariate file.\n"
"  --covar-name [...]      : Specify covariate(s) in --covar file by name.\n"
"                            Separate multiple names with spaces or commas, and\n"
"                            use dashes to designate ranges.\n"
"  --covar-number [...]    : Specify covariate(s) in --covar file by index.\n"
	       );
    help_print("within\tmwithin", &help_ctrl, 0,
"  --within [f] <keep-NA>  : Specify initial cluster assignments.\n"
"  --mwithin [n]           : Load cluster assignments from column n+2.\n"
	       );
    help_print("loop-assoc", &help_ctrl, 0,
"  --loop-assoc [f] <keep-NA>    : Run specified case/control association\n"
"                                  commands once for each cluster in the file,\n"
"                                  using cluster membership as the phenotype.\n"
	       );
#ifndef STABLE_BUILD
    help_print("set\tsubset\tset-collapse-all\tmake-set-collapse-all\tcomplement-sets\tmake-set-complement-all\tmake-set\tmake-set-border\tborder\tmake-set-collapse-group\t--make-set-complement-group", &help_ctrl, 0,
"  --set [filename]              : Load sets from a .set file.\n"
"  --subset [filename]           : Throw out sets not named in the given file.\n"
"  --set-collapse-all [set name] : Merge all sets.\n"
"  --complement-sets             : Invert all sets.  (Names gain 'C_' prefixes.)\n"
"  --make-set-complement-all [s] : --set-collapse-all + inversion.\n"
"  --make-set [filename]         : Define sets from a list of named bp ranges.\n"
"  --make-set-border [kbs]       : Stretch regions in --make-set file.\n"
"  --make-set-collapse-group     : Define sets from groups instead of sets in\n"
"                                  --make-set file.\n"
	       );
#endif
    help_print("keep\tremove\tkeep-fam\tremove-fam", &help_ctrl, 0,
"  --keep [fname]   : Exclude all individuals not named in the file.\n"
"  --remove [fname] : Exclude all individuals named in the file.\n"
"  --keep-fam [fn]  : Exclude all families not named in the file.\n"
"  --remove-fam [f] : Exclude all families named in the file.\n"
	       );
    help_print("extract\texclude", &help_ctrl, 0,
"  --extract [file] : Exclude all variants not named in the file.\n"
"  --exclude [file] : Exclude all variants named in the file.\n"
	       );
    help_print("keep-clusters\tkeep-cluster-names\tremove-clusters\tremove-cluster-names", &help_ctrl, 0,
"  --keep-clusters [filename]          : These can be used individually or in\n"
"  --keep-cluster-names [name(s)...]     combination to define a list of\n"
"                                        clusters to keep; all individuals not\n"
"                                        in a cluster in that list are then\n"
"                                        excluded.  Use spaces to separate\n"
"                                        cluster names for --keep-cluster-names.\n"
"  --remove-clusters [filename]        : Exclude all clusters named in the file.\n"
"  --remove-cluster-names [name(s)...] : Exclude the named clusters.\n"
	       );
#ifndef STABLE_BUILD
    help_print("set\tmake-set\tgene\tgene-all", &help_ctrl, 0,
"  --gene [sets...] : Exclude variants not in a set named on the command line.\n"
"  --gene-all       : Exclude variants which aren't a member of any set.  (PLINK\n"
"                     1.07 automatically did this under some circumstances.)\n"
	       );
#endif
    help_print("filter-attrib\tfilter-attrib-indiv", &help_ctrl, 0,
"  --filter-attrib [f] {att lst} : Given a file assigning attributes to\n"
"  --filter-attrib-indiv [f] {a}   variants, and a comma-delimited list (with no\n"
"                                  whitespace) of attribute names, remove\n"
"                                  variants/individuals which are either missing\n"
"                                  from the file or don't have any of the listed\n"
"                                  attributes.  If some attribute names in the\n"
"                                  list are preceded by '-', they are treated as\n"
"                                  \"negative match conditions\" instead: variants\n"
"                                  with all the negative match attributes are\n"
"                                  removed.\n"
"                                  The first character in the list cannot be a\n"
"                                  '-', due to how command-line parsing works;\n"
"                                  add a comma in front to get around this.\n"
	       );
    help_print("chr\tnot-chr\tchr-excl\tfrom-bp\tto-bp\tfrom-kb\tto-kb\tfrom-mb\tto-mb", &help_ctrl, 0,
"  --chr [chrs...]  : Exclude all variants not on the given chromosome(s).\n"
"                     Valid choices for humans are 0 (unplaced), 1-22, X, Y, XY,\n"
"                     and MT.  Separate multiple chromosomes with spaces and/or\n"
"                     commas, and use a dash (no adjacent spaces permitted) to\n"
"                     denote a range, e.g. '--chr 1-4, 22, xy'.\n"
"  --not-chr [...]  : Reverse of --chr (exclude variants on listed chromosomes).\n"
	       );
    help_print("autosome\tautosome-xy\tchr\tnot-chr\tchr-excl", &help_ctrl, 0,
"  --autosome       : Exclude all non-autosomal variants.\n"
"  --autosome-xy    : Exclude all non-autosomal variants, except those with\n"
"                     chromosome code XY (pseudo-autosomal region of X).\n"
	       );
    help_print("from\tto\tsnp\twindow\tfrom-bp\tto-bp\tfrom-kb\tto-kb\tfrom-mb\tto-mb\texclude-snp\textract-snp", &help_ctrl, 0,
"  --from [var ID]  : Use ID(s) to specify a variant range to load.  When used\n"
"  --to   [var ID]    together, both variants must be on the same chromosome.\n"
"  --snp  [var ID]  : Specify a single variant to load.\n"
"  --exclude-snp [] : Specify a single variant to exclude.\n"
"  --window  [kbs]  : With --snp or --exclude-snp, loads/excludes all variants\n"
"                     within half the specified kb distance of the named one.\n"
"  --from-bp [pos]  : Use physical position(s) to define a variant range to\n"
"  --to-bp   [pos]    load.  --from-kb/--to-kb/--from-mb/--to-mb allow decimal\n"
"    ...              values.  You must also specify a single chromosome (using\n"
"                     e.g. --chr) when using these flags.\n"
	       );
    help_print("snps\texclude-snps", &help_ctrl, 0,
"  --snps [var IDs...]  : Use IDs to specify variant range(s) to load or\n"
"  --exclude-snps [...]   exclude.  E.g. '--snps rs1111-rs2222, rs3333, rs4444'.\n"
	       );
    help_print("thin", &help_ctrl, 0,
"  --thin [p]       : Remove variants at random (p = retention probability).\n"
	       );
    help_print("bp-space", &help_ctrl, 0,
"  --bp-space [bps] : Remove variants so that each pair is no closer than the\n"
"                     given bp distance.\n"
	       );
    help_print("filter\tmfilter", &help_ctrl, 0,
"  --filter [f] [val(s)...] : Exclude all individuals without a 3rd column entry\n"
"                             in the given file matching one of the given\n"
"                             value(s).\n"
"  --mfilter [n]            : Match against (n+2)th column instead.\n"
	       );
    help_print("geno\tmind", &help_ctrl, 0,
"  --geno {val}     : Exclude variants with missing call frequencies greater\n"
"                     than a threshold (default 0.1).  (Note that the default\n"
"                     threshold is only applied if --geno is invoked without a\n"
"                     parameter; when --geno is not invoked, no per-variant\n"
"                     missing call frequency ceiling is enforced at all.  Other\n"
"                     inclusion/exclusion default thresholds work the same way.)\n"
"  --mind {val}     : Exclude individuals with missing call frequencies greater\n"
"                     than a threshold (default 0.1).\n"
	       );
    help_print("prune", &help_ctrl, 0,
"  --prune          : Remove individuals with missing phenotypes.\n"
	       );
    help_print("maf\tmax-maf", &help_ctrl, 0,
"  --maf {val}      : Exclude variants with minor allele frequency lower than a\n"
"                     threshold (default 0.01).\n"
"  --max-maf [val]  : Exclude variants with MAF greater than the threshold.\n"
	       );
    help_print("maf-succ", &help_ctrl, 0,
"  --maf-succ       : Rule of succession MAF estimation (used in EIGENSTRAT).\n"
"                     Given j observations of one allele and k >= j observations\n"
"                     of the other, infer a MAF of (j+1) / (j+k+2), rather than\n"
"                     the default j / (j+k).\n"
	       );
    help_print("read-freq\tupdate-freq", &help_ctrl, 0,
"  --read-freq [fn] : Estimate MAFs and heterozygote frequencies from the given\n"
"                     --freq[x] report, instead of the input fileset.\n"
	       );
    help_print("hwe", &help_ctrl, 0,
"  --hwe {val}      : Exclude variants with Hardy-Weinberg equilibrium exact\n"
"                     test p-values below a threshold (default 0.001).\n"
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
"  --make-founders <require-2-missing> <first> : Clear parental IDs for those\n"
"                                                with 1+ missing parent(s).\n"
	       );
    help_print("recode\trecode-allele", &help_ctrl, 0,
"  --recode-allele [f]  : With --recode A or --recode AD, count alleles named in\n"
"                         the file (instead of the minor allele).\n"
	       );
    help_print("set-hh-missing", &help_ctrl, 0,
"  --set-hh-missing : Cause --make-bed and --recode to set heterozygous haploid\n"
"                     genotypes to missing.\n"
	       );
    help_print("output-missing-genotype\toutput-missing-phenotype", &help_ctrl, 0,
"  --output-missing-genotype [ch] : Set the code used to represent missing\n"
"                                   genotypes in new filesets (normally the\n"
"                                   --missing-genotype value).\n"
"  --output-missing-phenotype [n] : Set the value used to represent missing\n"
"                                   phenotypes in new filesets (normally the\n"
"                                   --missing-phenotype value).\n"
	       );
    help_print("update-chr\tupdate-cm\tupdate-map\tupdate-name", &help_ctrl, 0,
"  --update-chr  [f] {chrcol} {IDcol}  {skip} : Update variant chromosome codes.\n"
"  --update-cm   [f] {cmcol}  {IDcol}  {skip} : Update centimorgan positions.\n"
"  --update-map  [f] {bpcol}  {IDcol}  {skip} : Update variant bp positions.\n"
"  --update-name [f] {newcol} {oldcol} {skip} : Update variant IDs.\n"
	       );
    help_print("update-alleles", &help_ctrl, 0,
"  --update-alleles [f]     : Update variant allele codes.\n"
	       );
    help_print("allele1234\talleleACGT\talleleacgt", &help_ctrl, 0,
"  --allele1234 <multichar> : Interpret/recode A/C/G/T alleles as 1/2/3/4.\n"
"                             With 'multichar', converts all A/C/G/Ts in allele\n"
"                             names to 1/2/3/4s.\n"
"  --alleleACGT <multichar> : Reverse of --allele1234.\n"
	       );
    help_print("update-ids\tupdate-parents\tupdate-sex", &help_ctrl, 0,
"  --update-ids [f]     : Update individual IDs.\n"
"  --update-parents [f] : Update parental IDs.\n"
"  --update-sex [f]     : Update individual sexes.\n"
	       );
    help_print("flip\tflip-subset", &help_ctrl, 0,
"  --flip [filename]    : Flip alleles (A<->T, C<->G) for SNP IDs in the file.\n"
#ifndef STABLE_BUILD
"  --flip-subset [fn]   : Only apply --flip to indivs in the --flip-subset file.\n"
#endif
	       );
    help_print("keep-allele-order\tmake-bed\tmerge\tbmerge\tmerge-list", &help_ctrl, 0,
"  --keep-allele-order  : Keep the allele order defined in the .bim file,\n"
"                         instead of forcing A2 to be the major allele.\n"
	       );
    help_print("a1-allele\treference-allele\tupdate-ref-allele\ta2-allele", &help_ctrl, 0,
"  --a1-allele [f] {a1col} {IDcol} {skip} : Force alleles in the file to A1.\n"
"  --a2-allele [f] {a2col} {IDcol} {skip} : Force alleles in the file to A2.\n"
	       );
    help_print("indiv-sort\tmerge\tbmerge\tmerge-list", &help_ctrl, 0,
"  --indiv-sort [m] : Specify family/individual ID sort order.  The following\n"
"                     three modes are currently supported:\n"
"                     * 'none'/'0' keeps individuals in the order they were\n"
"                       loaded.  This is the default for non-merge operations.\n"
"                     * 'natural'/'n' invokes \"natural sort\", e.g. 'id2' <\n"
"                       'ID3' < 'id10'.  This is the default when merging.\n"
"                     * 'ascii'/'a' sorts in ASCII order, e.g. 'ID3' < 'id10' <\n"
"                       'id2'.\n"
"                     For now, only --make-bed and --merge/--bmerge/--merge-list\n"
"                     respect this flag.\n"
	       );
    help_print("with-phenotype\tdummy-coding\twrite-covar", &help_ctrl, 0,
"  --with-phenotype <no-parents> <no-sex | female-2> : Include more individual\n"
"                                                      info in new .cov file.\n"
"  --dummy-coding {N} <no-round> : Split categorical variables (n categories,\n"
"                                  2 < n <= N, default N is 49) into n-1 binary\n"
"                                  dummy variables when writing covariate file.\n"
	       );
    help_print("merge\tbmerge\tmerge-list\tmerge-mode", &help_ctrl, 0,
"  --merge-mode [n]   : Adjust --merge/--bmerge/--merge-list behavior based on a\n"
"                       numeric code.\n"
"                       1 (default) = difference -> missing\n"
"                       2 = only overwrite originally missing calls\n"
"                       3 = only overwrite when nonmissing in new file\n"
"                       4/5 = never overwrite and always overwrite, respectively\n"
"                       6 = report all mismatching calls without merging\n"
"                       7 = report mismatching nonmissing calls without merging\n"
	       );
    help_print("merge\tbmerge\tmerge-list\tmerge-mode\tmerge-equal-pos", &help_ctrl, 0,
"  --merge-equal-pos  : Merge variants with different names but identical\n"
"                       positions.\n"
	       );
#ifndef STABLE_BUILD
    help_print("r\tr2\tld-window-r2\tld-window\tld-window-kb\tld-snp\tld-snps\tld-snp-list", &help_ctrl, 0,
"  --ld-window [ct+1] : Set --r/--r2 max site ct pairwise distance (usually 10).\n"
"  --ld-window-kb [x] : Set --r/--r2 max kb pairwise distance (usually 200).\n"
"  --ld-window-r2 [x] : Set threshold for --r2 report inclusion (usually 0.2).\n"
"  --ld-snp [var ID]  : Set first variant in all --r/--r2 pairs.\n"
"  --ld-snps [vID...] : Restrict first --r/--r2 variant to the given ranges.\n"
"  --ld-snp-list [f]  : Restrict first --r/--r2 var. to those named in the file.\n"
	       );
#endif
    help_print("indep\tindep-pairwise\tld-xchr", &help_ctrl, 0,
"  --ld-xchr [code]   : Specify X chromosome model for --indep[-pairwise].\n"
"                       1 (default) = males coded 0/1, females 0/1/2 (A1 dosage)\n"
"                       2 = males coded 0/2\n"
"                       3 = males coded 0/2, but females given double weighting\n"
	       );
    help_print("distance-exp\texponent\tdistance", &help_ctrl, 0,
"  --distance-exp [x] : When computing genomic distances, assign each variant a\n"
"                       weight of (2q(1-q))^{-x}, where q is the inferred MAF.\n"
"                       (Use --read-freq if you want to explicitly specify some\n"
"                       or all of the MAFs.)\n"
	       );
    help_print("read-dists\tload-dists\tibs-test\tgroupdist\tregress-distance\tcluster\tneighbour\tneighbor", &help_ctrl, 0,
"  --read-dists [dist file] {id file} : Load a triangular binary distance matrix\n"
"                                       instead of recalculating from scratch.\n"
	       );
    help_print("ppc-gap\tmin\tmax\tgenome\tZ-genome", &help_ctrl, 0,
"  --ppc-gap [val]    : Minimum number of base pairs, in thousands, between\n"
"                       informative pairs of markers used in --genome PPC test.\n"
"                       500 if unspecified.\n"
"  --min [cutoff]     : Specify minimum PI_HAT for inclusion in --genome report.\n"
"  --max [cutoff]     : Specify maximum PI_HAT for inclusion in --genome report.\n"
	       );
    help_print("homozyg\thomozyg-match\tpool-size", &help_ctrl, 0,
"  --homozyg-match [] : Set minimum concordance across jointly homozygous sites\n"
"                       for a pairwise allelic match to be declared.\n"
"  --pool-size [ct]   : Set minimum size of pools in '--homozyg group' report.\n"
	       );
    help_print("read-genome\tcluster\tneighbour\tneighbor", &help_ctrl, 0,
"  --read-genome [fn] : Load --genome report for --cluster/--neighbour, instead\n"
"                       of recalculating IBS and PPC test p-values from scratch.\n"
	       );
    help_print("ppc\tmc\tmcc\tK\tk\tibm\tcluster", &help_ctrl, 0,
"  --ppc [p-val]    : Specify minimum PPC test p-value within a cluster.\n"
"  --mc [max size]  : Specify maximum cluster size.\n"
"  --mcc [c1] [c2]  : Specify maximum case and control counts per cluster.\n"
"  --K [min count]  : Specify minimum cluster count.\n"
"  --ibm [val]      : Specify minimum identity-by-missingness.\n"
	       );
    help_print("match\tmatch-type\tqmatch\tqt\tcluster", &help_ctrl, 0,
"  --match [f] {mv} : Use covariate values to restrict clustering.  Without\n"
"                     --match-type, two individuals can only be in the same\n"
"                     cluster if all covariates match.  The optional second\n"
"                     parameter specifies a covariate value to treat as missing.\n"
"  --match-type [f] : Refine interpretation of --match file.  The --match-type\n"
"                     file is expected to be a single line with as many entries\n"
"                     as the --match file has covariates; '0' entries specify\n"
"                     'negative matches' (i.e. individuals with equal covariate\n"
"                     values cannot be in the same cluster), '1' entries specify\n"
"                     'positive matches' (default), and '-1' causes the\n"
"                     corresponding covariate to be ignored.\n"
"  --qmatch [f] {m} : Force all members of a cluster to have similar\n"
"  --qt [fname]       quantitative covariate values.  The --qmatch file contains\n"
"                     the covariate values, while the --qt file is a list of\n"
"                     nonnegative tolerances (and '-1's marking covariates to\n"
"                     skip).\n"
	       );
#ifndef NOLAPACK
    help_print("cluster\tmds-plot\tmds-cluster", &help_ctrl, 0,
"  --mds-plot [dims] <by-cluster> <eigvals> : Multidimensional scaling analysis.\n"
"                                             Requires --cluster.\n"
	       );
#endif
    help_print("cell\tmodel", &help_ctrl, 0,
"  --cell [thresh]  : Skip some --model tests when a contingency table entry is\n"
"                     smaller than the given threshold.\n"
	       );
    help_print("linear\tlogistic\tcondition\tcondition-list\tparameters\ttests\ttest-all\tvif\txchr-model", &help_ctrl, 0,
"  --condition [var ID] <dominant | recessive> : Add one variant as a --linear\n"
"                                                or --logistic covariate.\n"
"  --condition-list [f] <dominant | recessive> : Add variants named in the file\n"
"                                                as --linear/--logistic covs.\n"
"  --parameters [...]  : Include only the given covariates/interactions in the\n"
"                        --linear/--logistic models, identified by a list of\n"
"                        1-based indices and/or ranges of them.\n"
"  --tests <all> {...} : Perform a (joint) test on the specified term(s) in the\n"
"                        --linear/--logistic model, identified by 1-based\n"
"                        indices and/or ranges of them.  If permutation was\n"
"                        requested, it is based on this test.\n"
"                        * Note that, when --parameters is also present, the\n"
"                          indices refer to the terms remaining AFTER pruning by\n"
"                          --parameters.\n"
"                        * You can use '--tests all' to include all terms.\n"
"  --vif [max VIF]     : Set VIF threshold for --linear/--logistic\n"
"                        multicollinearity check (default 50).\n"
"  --xchr-model [code] : Set the X chromosome --linear/--logistic model.\n"
"                        0 = skip sex and haploid chromosomes\n"
"                        1 (default) = add sex as a covariate on X chromosome\n"
"                        2 = code male genotypes 0/2 instead of 0/1\n"
"                        3 = test for interaction between genotype and sex\n"
	       );
    help_print("lasso\tlasso-select-covars", &help_ctrl, 0,
"  --lasso-select-covars {cov(s)...} : Subject some or all covariates to LASSO\n"
"                                      model selection.\n"
	       );
    help_print("adjust\tgc\tlog10\tqq-plot", &help_ctrl, 0,
"  --adjust <gc> <log10> <qq-plot>   : Report some multiple-testing corrections.\n"
	       );
    help_print("adjust\tlambda", &help_ctrl, 0,
"  --lambda [val]   : Set genomic control lambda for --adjust.\n"
	       );
    help_print("ci", &help_ctrl, 0,
"  --ci [size]      : Report confidence intervals for odds ratios.\n"
	       );
    help_print("pfilter", &help_ctrl, 0,
"  --pfilter [val]  : Filter out association test results with higher p-values.\n"
	       );
    help_print("aperm", &help_ctrl, 0,
"  --aperm [min perms - 1] {max perms} {alpha} {beta} {init interval} {slope} :\n"
"    Set up to six parameters controlling adaptive permutation tests.\n"
"    * The first two control the minimum and maximum number of permutations that\n"
"      may be run for each variant; default values are 5 and 1000000.\n"
"    * The next two control the early termination condition.  A\n"
"      100% * (1 - beta/2T) confidence interval is calculated for each empirical\n"
"      p-value, where T is the total number of variants; whenever this\n"
"      confidence interval doesn't contain alpha, the variant is exempted from\n"
"      further permutation testing.  Default values are 0 and 0.0001.\n"
"    * The last two control when the early termination condition is checked.  If\n"
"      a check occurs at permutation #p, the next check occurs after\n"
"      [slope]p + [init interval] more permutations (rounded down).  Default\n"
"      initial interval is 1, and default slope is 0.001.\n"
	       );
    help_print("mperm-save\tmperm-save-all", &help_ctrl, 0,
"  --mperm-save     : Save best max(T) permutation test statistics.\n"
"  --mperm-save-all : Save all max(T) permutation test statistics.\n"
	       );
    help_print("fast-epistasis\tepistasis\tgap\tepi1\tepi2", &help_ctrl, 0,
"  --gap [kbs]      : Set '--fast-epistasis case-only' min. gap (default 1000).\n"
"  --epi1 [p-value] : Set --[fast-]epistasis reporting threshold (def. 0.0001).\n"
"  --epi2 [p-value] : Set threshold for contributing to SIG_E count (def. 0.01).\n"
	       );
    help_print("parallel\tgenome-lists", &help_ctrl, 0,
"  --parallel [k] [n] : Divide the output matrix into n pieces, and only compute\n"
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
"  --memory [val]     : Set size, in MB, of initial workspace malloc attempt.\n"
	       );
    help_print("threads\tthread-num\tnum_threads", &help_ctrl, 0,
"  --threads [val]    : Set maximum number of concurrent threads.\n"
	       );
    help_print("d\tsnps", &help_ctrl, 0,
"  --d [char]         : Change variant/covariate range delimiter (normally '-').\n"
	       );
    help_print("seed", &help_ctrl, 0,
"  --seed [val...]    : Set random number seed(s).  Each value must be an\n"
"                       integer between 0 and 4294967295 inclusive.\n"
	       );
    help_print("perm-batch-size", &help_ctrl, 0,
"  --perm-batch-size [val] : Set number of permutations per batch in QT\n"
"                            permutation tests.\n"
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
    help_print("cnv-exclude-off-by-1", &help_ctrl, 0,
"  --cnv-exclude-off-by-1   : Exclude .cnv segments where the terminal .cnv.map\n"
"                             entry is off by 1.\n"
	       );
    help_print("cnv-test-window\tcnv-test", &help_ctrl, 0,
"  --cnv-test-window [size] : Specify window size (in kb) for CNV assoc. test.\n"
	       );
    help_print("cnv-count\tcnv-indiv-perm\tcnv-enrichment-test", &help_ctrl, 0,
"  --cnv-count [fname]      : Specify region list for --cnv-indiv-perm\n"
"                             (optional) or --cnv-enrichment-test (required).\n"
	       );
#endif
    if (!param_ct) {
      fputs(
"\nFor further documentation and support, consult the main webpage\n"
#ifdef PLINK_BUILD
"([TBD]) and/or the mailing list ([TBD]).\n"
#else
"(https://www.cog-genomics.org/wdist ) and/or the wdist-users mailing list\n"
"(https://groups.google.com/d/forum/wdist-users ).\n"
#endif
, stdout);
    }
  } while (help_ctrl.iters_left--);
  if (help_ctrl.unmatched_ct) {
    net_unmatched_ct = help_ctrl.unmatched_ct;
    printf("\nNo help entr%s for", (help_ctrl.unmatched_ct == 1)? "y" : "ies");
    col_num = (help_ctrl.unmatched_ct == 1)? 17 : 19;
    arg_uidx = 0;
    while (help_ctrl.unmatched_ct) {
      arg_uidx = next_unset_unsafe(help_ctrl.all_match_arr, arg_uidx);
      help_ctrl.unmatched_ct--;
      if (help_ctrl.unmatched_ct) {
	if (net_unmatched_ct == 2) {
	  if (help_ctrl.param_lens[arg_uidx] + col_num > 76) {
	    putchar('\n');
	    col_num = 2 + help_ctrl.param_lens[arg_uidx];
	  } else {
	    putchar(' ');
	    col_num += 3 + help_ctrl.param_lens[arg_uidx];
	  }
	  putchar('\'');
	  fputs(argv[arg_uidx], stdout);
	  putchar('\'');
	} else {
	  if (help_ctrl.param_lens[arg_uidx] + col_num > 75) {
	    putchar('\n');
	    col_num = 3 + help_ctrl.param_lens[arg_uidx];
	  } else {
	    putchar(' ');
	    col_num += 4 + help_ctrl.param_lens[arg_uidx];
	  }
	  putchar('\'');
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
	putchar((help_ctrl.param_lens[arg_uidx] + col_num > 75)? '\n' : ' ');
	putchar('\'');
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

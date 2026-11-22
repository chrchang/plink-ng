

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "plink.h"
#include "helper.h"
#include "options.h"
#include "perm.h"
#include "stats.h"
#include "linear.h"
#include "logistic.h"


extern Plink * PP;
class RCount
{
public:

  RCount(Plink * p_, map<Range,int2> * rl_)
  {
    P = p_;

    // not used now:
    rangeLookup = rl_;

    rval.resize(P->n,0);  // Potentially weighted score
    aval.resize(P->n,0);  // 0/1/2 for present rare allele
    gval.resize(P->n,0);  // 0/1 for genotyped or not
    nsnps = nalleles = 0;
    npc = 0;
    pcMode = false;
    domModel = true;
  }

  Plink * P;
  map<Range,int2> * rangeLookup;

  vector_t rval;
  vector<int> gval;
  vector<int> aval;

  // Current SNPs in window -> SNP specific counts
  map<int,vector_t > rwin;
  map<int,vector<int> > gwin;
  map<int,vector<int> > awin;

  double acnt, ucnt;
  int acnt2, ucnt2;
  int nsnps;
  int nalleles;
  int npc;

  bool pcMode;
  bool domModel;

  bool addSNP(int l);
  bool removeSNP(int l);
  bool setWindow(int chr, int bp);
  void displayWindow();
  void loadCovariate();
  void loadPCACovariate();
  void mainStats();
  bool ignorePosition();
};

bool RCount::ignorePosition()
{
  if ( nsnps < 1 )
    return true;
  return nalleles < 5;
}

void RCount::loadCovariate()
{
  // Place as last covariate
  for (int i=0; i<P->n; i++)
    {
      P->sample[i]->clist[ par::clist_number - 1 ] =
	gval[i] > 0 ?
	rval[i] / (double)(gval[i]) :
	0;
    }
}

void RCount::loadPCACovariate()
{

  // In place of simply counting all the rare SNPs, perform a PCA on the
  // rare SNP data matrix, then sum the standardized PCs above a certian
  // threshold (i.e. this way giving equal weight to equal independently
  // detected component of rare variation)

  vector<int> snplist;
  map<int,vector_t>::iterator i = rwin.begin();
  while ( i != rwin.end() )
    {
      snplist.push_back( i->first );
      ++i;
    }

  boolmatrix_t mask;
  matrix_t g;

  bool dominantModel = true;
  geno2matrix( snplist , g , mask , dominantModel );

  vector_t pc;
  matrix_t ps;
  matrix_t pv;
  matrix_t g0 = g;


  // Setting last flag to false implies no mean-centering
  // This version of PCA return U.W*.V' in 'g', where
  // W* is an editted eigen-value set, such that they
  // equal eithe 0 or 1

  bool meanCentre = par::elf_pcmode_2sided;

  int pcn = pca( g , mask , pc , ps , pv, meanCentre);
  int ntot = g[0].size();


  if ( ! par::elf_pcmode_2sided )
    {

      // Calculate score which is sum of squares of U.W*.V'
      vector_t sc(P->n,0);

      for (int i=0; i<P->n; i++)
	{
	  for (int p = 0 ; p < ntot ; p++ )
	    sc[i] += g[i][p] * g[i][p];
	}

      // Standardize and threshold at 4SD
      double m = 0;
      double ssq = 0;
      double v = 0;
      for (int i=0; i<P->n; i++)
	{
	  m += sc[i];
	  ssq += sc[i] * sc[i];
	}
      m /= P->n;
      ssq /= P->n;
      v = ssq - m*m;
      double sd = sqrt(v);

      // Now assign...
      for (int i=0; i<P->n; i++)
	{
	  double z = (sc[i]-m)/sd;
	  if ( z > 4) z = 4;
	  P->sample[i]->clist[ par::clist_number - 1 ] = z;
	}
    }
  else
    {

      // Resize covariate load
      par::clist_number -= npc;

      // Place as set of covariates, 1 past last covariate
      int start = par::clist_number;

      par::clist_number += pcn;

      // End
      int end = par::clist_number-1;

      P->clistname.resize( par::clist_number );

//       cout << "pcn, npc = " << pcn << " " << npc << "\n";
//       cout << "Loading covars : " << start << " to " << end << "\n";

      for (int i=0; i<P->n; i++)
	{
	  P->sample[i]->clist.resize( par::clist_number );

	  int k=0;
	  for (int j=start; j<=end; j++)
	    {
	      P->sample[i]->clist[ j ] = ps[i][k++];
	    }
	}

      // Track number of PCs so they can be accounted for in the next analysis
      npc = pcn;
    }

}

bool RCount::setWindow(int chr, int bp)
{

  // Find all SNPs with x kb of bp, and add to rwin, if not already in
  // there Also, keep track of what we have added, and remove any SNPs
  // that should no longer be in the window


  // Return true if window actually changes since last position
  bool changed = false;

  set<int> nwin;

  // Use a Range to lookup the SNPs in this range

  Range r;
  r.chr = chr;
  r.start = bp - (int)par::rarer_dist_threshold;
  r.stop = bp + (int)par::rarer_dist_threshold;

  // Start and stop sites for this range:
  int2 l2 = mapSNPs2Range( *PP , &r );


  /////////////////////////////
  // Need to add any new SNPs?

  if ( l2.p1 != -1 )
    for (int l = l2.p1 ; l <= l2.p2 ; l++ )
      {
	if ( addSNP(l) )
	  changed = true;
	nwin.insert(l);
      }


  ////////////////////////
  // Need to remove any?

  map<int, vector_t >::iterator iter = rwin.begin();
  set<int> toRemove;
  while ( iter != rwin.end() )
    {
      if ( nwin.find( iter->first ) == nwin.end() )
	toRemove.insert( iter->first );
      ++iter;
    }

  set<int>::iterator i2 = toRemove.begin();
  while( i2 != toRemove.end() )
    {
      if ( removeSNP( *i2 ) )
	changed = true;
      ++i2;
    }

  return changed;
}


void RCount::displayWindow()
{
  int rmin = 9999999;
  int rmax = -1;

  map<int, vector_t >::iterator iter = rwin.begin();
  while ( iter != rwin.end() )
    {
      if ( iter->first < rmin )
	rmin = iter->first;
      if ( iter->first > rmax )
	rmax = iter->first;
      ++iter;
    }

  acnt = 0;
  ucnt = 0;
  for (int i=0; i<P->n; i++)
    if ( P->sample[i]->pperson->aff )
      acnt += rval[i];
    else
      ucnt += rval[i];

  cout << "Window from "
       << P->locus[rmin]->name << "("
       << P->locus[rmin]->bp << ") to "
       << P->locus[rmax]->name << "("
       << P->locus[rmax]->bp << ") to "
       << rwin.size() << " SNPs with vals (A/U) "
       << acnt << " and "
       << ucnt
       << "\n";

}


void RCount::mainStats()
{

  nsnps = rwin.size();
  acnt = 0;
  ucnt = 0;
  acnt2 = 0;
  ucnt2 = 0;

  int acnt3 = 0;
  int ucnt3 = 0;

  for (int i=0; i<P->n; i++)
    {
      if ( ! P->sample[i]->missing )
	{
	  if ( par::bt )
	    {
	      if ( P->sample[i]->pperson->aff )
		{
		  acnt += rval[i];
		  ++acnt2;
		  acnt3 += aval[i];
		}
	      else
		{
		  ucnt += rval[i];
		  ++ucnt2;
		  ucnt3 += aval[i];
		}
	    }
	  else
	    {
	      ucnt += rval[i];
	      ++ucnt2;
	      ucnt3 += aval[i];
	    }
	}
    }

  // Number of low-frequency alleles in
  // this window

  nalleles = (int)acnt3 + (int)ucnt3;

  // The proportion of alleles that are LF

  if ( par::bt && acnt2>0)
    acnt /= (double)acnt2;
  if ( ucnt2>0 )
    ucnt /= (double)ucnt2;

}


bool RCount::addSNP(int l)
{

  if ( rwin.find(l) != rwin.end() )
    return false;


  if ( P->locus[l]->freq > par::rarer_maf_threshold )
    {
      return false;
    }

  vector_t r(P->n,0);
  vector<int> g(P->n,0);
  vector<int> a(P->n,0);

  double wt;
  if ( par::rare_test_weight1 )
    wt = 1/P->locus[l]->freq;

  CSNP * snp = P->SNP[l];

  for (int i=0; i<P->n; i++)
    {

      //////////////////////////////////////////
      // Get and parse genotypes

      bool one = snp->one[i];
      bool two = snp->two[i];

      // Skip if missing

      if ( one && !two )
	continue;

      if ( domModel )
	{
	  // Dominant coding

	  if ( ( ! one ) || ( ! two ) )
	    {
	      r[i] += par::rare_test_weight1 ? wt : 1 ;
	      ++a[i];
	    }

	  g[i] += 1;

	}
      else
	{

	  // Additive coding

	  if ( ! one )
	    {
	      r[i] += par::rare_test_weight1 ? wt : 1 ;
	      ++a[i];
	    }

	  if ( ! two )
	    {
	      r[i] += par::rare_test_weight1 ? wt : 1 ;
	      ++a[i];
	    }

	  g[i] += 2;
	}

      // Add to current total per person
      rval[i] += r[i];
      aval[i] += a[i];
      if ( domModel )
	gval[i] += 1;
      else
	gval[i] += 2;
    }

  rwin.insert(make_pair( l , r ));
  gwin.insert(make_pair( l , g ));
  awin.insert(make_pair( l , a ));
  return true;
}

bool RCount::removeSNP(int l)
{

  map<int, vector_t >::iterator iter = rwin.find(l);

  // If SNP never added to window, nothing to do

  if ( iter == rwin.end() )
    return false;

  map<int, vector<int> >::iterator giter = gwin.find(l);
  map<int, vector<int> >::iterator aiter = awin.find(l);

  for (int i=0; i<P->n; i++)
    {
      rval[i] -= iter->second[i];
      gval[i] -= giter->second[i];
      aval[i] -= aiter->second[i];
    }

  rwin.erase(iter);
  gwin.erase(giter);
  awin.erase(aiter);
  return true;
}


// Other output helper functions
void displayScoresPerson(ofstream & O, RCount & rc)
{
  for (int i = 0 ; i < PP->n ; i++ )
    {
      O << setw(par::pp_maxfid ) << PP->sample[i]->fid << " "
	<< setw(par::pp_maxiid ) << PP->sample[i]->iid << " ";
      if ( PP->sample[i]->missing )
	O << "NA" << "\t" << "NA" << "\t" << "NA" << "\n";
      else
	O << PP->sample[i]->phenotype << "\t"
	  << PP->sample[i]->clist[ par::clist_number - 1 ] << "\t"
	  << rc.aval[i] << "\t"
	  << rc.gval[i] << "\n";
    }
}

void displayScoresRegion(ofstream & O, RCount & rc)
{

  map<int,vector<int> >::iterator i = rc.awin.begin();
  while ( i != rc.awin.end() )
    {
      int count = 0;
      for ( int k = 0 ; k < i->second.size(); k++)
	count += i->second[k];

      O << setw(4) << PP->locus[ i->first ]->chr << " "
	<< setw(par::pp_maxsnp ) << PP->locus[ i->first ]->name << " "
	<< setw(12) << PP->locus[ i->first ]->bp << " "
	<< setw(12) << PP->locus[ i->first ]->freq << " "
	<< setw(12) << PP->locus[ i->first ]->allele1 << " "
	<< setw(12) << count << "\n";
      ++i;
    }
}


void Plink::permTestRareDistribution(Perm & perm)
{

  printLOG("Testing for Enrichment of Low Frequency variants ");

  if ( par::rare_test_weight1 )
    printLOG(" (1/MAF weighting, ");
  else
    printLOG(" ( MAF < "
	     +dbl2str(par::rarer_maf_threshold)
	     +", ");

  printLOG("within "
	   +int2str(int(par::rarer_dist_threshold/1000))
	   +" kb)\n");



  ////////////////////////////
  // Use last covariate slot

  par::assoc_glm_without_main_snp = true;

  par::clist = true;

  if ( !par::elf_pcmode_2sided )
    {
      ++par::clist_number;
      clistname.push_back("RCNT");
      for (int i=0; i<n; i++)
	sample[i]->clist.push_back(0);
    }


  // NOTE: Not used now
  map<Range,int2> ranges;


  ///////////////////////
  // Original

  vector_t original = testRareDistribution(perm, true, ranges);


  ///////////////////////
  // Set up permutations

  perm.setTests( original.size() );
  perm.setPermClusters(*this);
  perm.originalOrder();

  if ( ! par::permute )
    return;

  if (par::mperm_rank)
    perm.setOriginalRanking(original);


  //////////////////////
  // Begin permutations

  bool finished = false;
  while(!finished)
    {
      perm.permuteInCluster();
      vector_t pr = testRareDistribution(perm,false, ranges);
      finished = perm.update(pr,original);
    }

  if (!par::silent)
    cout << "\n\n";



  /////////////////////////////////////////////////////////////////////
  // Write results to file

  ofstream ASC;
  string f;

  if (par::adaptive_perm) f = par::output_file_name + ".elf.perm";
  else f = par::output_file_name + ".elf.mperm";

  ASC.open(f.c_str(),ios::out);
  ASC.precision(4);

  printLOG("Writing permutation association results to [ " + f + " ] \n");

  ASC << setw(4) << "CHR" << " "
      << setw(par::pp_maxsnp)<< "SNP" << " "
      << setw(12)<< "STAT" << " "
      << setw(12) << "EMP1" << " ";
  if (par::adaptive_perm)
    ASC << setw(12)<< "NP" << " ";
  else if ( par::mperm_rank )
    ASC << setw(12)<< "EMP3" << " "
        << setw(12)<< "RANK" << " ";
  else
    ASC << setw(12)<< "EMP2" << " ";
  ASC << "\n";


  for (int l=0; l< original.size(); l++)
    {

      // Skip?, if filtering p-values
      if ( par::pfilter && perm.pvalue(l) > par::pfvalue )
        continue;

//       ASC << setw(4)  << locus[l]->chr << " "
//           << setw(par::pp_maxsnp) << locus[l]->name << " ";


      ASC << setw(8) << l << " ";

      ASC << setw(12) << original[l]  << " "
          << setw(12) << perm.pvalue(l) << " ";

      if (par::adaptive_perm)
        ASC << setw(12) << perm.reps_done(l) << " ";
      else if ( par::mperm_rank )
        ASC << setw(12) << perm.max_pvalue(l) << " "
            << setw(12) << perm.rank(l) << " ";
      else
        ASC << setw(12) << perm.max_pvalue(l) << " ";

      ASC << "\n";

    }

  ASC.close();

}



vector_t Plink::testRareDistribution(Perm & perm , bool disp,
				     map<Range,int2> & ranges)
{

  /////////////////////////////////////////////////////////////////////
  // Write results to file

  ofstream OUT;
  OUT.precision(4);

  ofstream SUM;
  const double pthresh = 0.01;
  bool one_sided = true;

  ofstream SDET_SNP;
  ofstream SDET_IND;

  if ( disp )
    {

      string f = par::output_file_name + ".elf";
      OUT.open(f.c_str(),ios::out);
      printLOG("Writing results to [ " + f + " ]\n");

      OUT  << setw(4) << "CHR" << " "
	   << setw(12) << "BP1" << " "
	   << setw(12) << "BP2" << " "
	   << setw(12) << "BP" << " "
	   << setw(6) << "NSNP" << " "
	   << setw(8) << "NALLELE" << " ";
      if ( par::bt )
	{OUT << setw(8) << "ACNT" << " "
	     << setw(8) << "UCNT" << " ";
	if ( ! par::elf_pcmode_2sided )
	  OUT << setw(10) << "OR" << " ";
	}
      else
	{
	  OUT << setw(8) << "CNT" << " ";
	  if ( ! par::elf_pcmode_2sided )
	    OUT << setw(10) << "BETA" << " ";
	}
      OUT << setw(10) << "CHISQ" << " "
	  << setw(10) << "P" << "\n";

    }


  // Use regression model: put # of rare variants per individual as a
  // covariate, and use glmAssoc()


  // We do not know how many results we will obtain to start off with

  vector_t results;


  RCount rc(this,&ranges);
  if ( par::elf_pcmode )
    rc.pcMode = true;

  vector_t b;
  double chisq;
  double pvalue;
  int srange_cnt = 0;
  bool inRange = false;

  int startChromosome = locus[ 0 ]->chr;
  int finalChromosome = locus[ nl_all - 1]->chr;

  for (int chr = startChromosome ; chr <= finalChromosome; chr++)
    {

      int bpstart = scaffold[chr].bpstart;
      int bpstop = scaffold[chr].bpstop;

      for ( int bp = bpstart; bp <= bpstop; bp += par::rarer_interval )
	{

	  bool windowMoved = rc.setWindow(chr,bp);

	  //rc.displayWindow();

	  // If no new SNPs have been added or removed from window,
	  // then just serve up the same results as last time

	  if ( ! windowMoved )
	    {
	      continue;
	    }


	  rc.mainStats();


	  // Enough to be bothering with?

	  if ( rc.ignorePosition() )
	    {
	      continue;
	    }

	  // Perform actual test

	  if ( rc.pcMode )
	    rc.loadPCACovariate();
	  else
	    rc.loadCovariate();

	  glmAssoc( false , perm );

	  // Get results

	  double beta;


	  if ( ! par::elf_pcmode_2sided )
	    {

	      model->testParameter = par::clist_number;
	      b = model->getCoefs();
	      chisq = model->getStatistic();
	      pvalue = chiprobP(chisq,1);
	      beta = b[ par::clist_number ];

	    }
	  else
	    {
	      b = model->getCoefs();
	      beta = 1;  // no direction
	      vector_t h; // dim = number of fixes (to =0)
	      matrix_t H; // row = number of fixes; cols = np
	      h.resize(rc.npc,0);
	      sizeMatrix(H,rc.npc,model->getNP());

	      int startpc = par::clist_number - rc.npc;

	      for (int i=0; i<rc.npc; i++)
		H[i][startpc + i + 1 ] = 1;

	      chisq = model->isValid() ? model->linearHypothesis(H,h) : 0;
	      pvalue = model->isValid() ? chiprobP(chisq, rc.npc) : -1;
	    }


	  // Permutation test is 1-sided

	  if ( beta < 0 )
	    results.push_back( 0 );
	  else
	    results.push_back( chisq );


	  // Clean up

	  delete model;


	  // Write results to a file?

	  if ( disp )
	    {

	      double coef = par::bt ? exp( beta ) : beta ;

	      int bp1 = bp - (int)par::rarer_dist_threshold < bpstart ?
		bpstart : bp - (int)par::rarer_dist_threshold;
	      int bp2 = bp + (int)par::rarer_dist_threshold > bpstop ?
		bpstop   : bp + (int)par::rarer_dist_threshold;

	      OUT  << setw(4) << chr << " "
		   << setw(12) << bp1 << " "
		   << setw(12) << bp2 << " "
		   << setw(12) << (int)((bp1+bp2)/2.0) << " "
		   << setw(6) << rc.nsnps << " "
		   << setw(8) << rc.nalleles << " ";
	      if ( par::bt ) OUT << setw(8) << rc.acnt << " ";
	      OUT << setw(8) << rc.ucnt << " ";

	      if ( ! par::elf_pcmode_2sided )
		OUT << setw(10) << coef << " ";
	      OUT  << setw(10) << chisq << " "
		  << setw(10) << pvalue << "\n";

	      OUT.flush();

	      if ( par::rare_test_print_details &&
		   int2str(chr)+":"+int2str( (int)((bp1+bp2)/2.0)) == par::rare_test_print_details_snp )
		{

		  printLOG("Printing details for region around " +  par::rare_test_print_details_snp + "\n");
		  SDET_SNP.open( ( par::output_file_name+".elf.det."
				   + par::rare_test_print_details_snp + ".snp" ).c_str() , ios::out );
		  SDET_IND.open( ( par::output_file_name+".elf.det."
				   + par::rare_test_print_details_snp + ".ind" ).c_str() , ios::out );

		  //////////////////////////////////
		  // Print scores per person, and per SNP

		  displayScoresPerson( SDET_IND , rc );
		  SDET_IND.close();
		  displayScoresRegion( SDET_SNP , rc );
		  SDET_SNP.close();

		}

	    } // end if verbose display mode


	} // Next window location
    }


  // Finished, close any open streams

  if ( disp )
    {
      OUT.close();
    }

  return results;

}



void Plink::displayRareRange()
{

  map<string, set<Range> > ranges = readRange( par::rare_test_score_range_file );

  printLOG("Reading ELF results file from [ "
	   + par::rare_test_score_results_file + " ]\n");

  printLOG("Reading ELF ranges from [ "
	   + par::rare_test_score_range_file + " ]\n");

  checkFileExists( par::rare_test_score_results_file );

  ifstream IN;
  IN.open( par::rare_test_score_results_file.c_str() );

  // Read first row

  int pcol = -1;
  int bpcol = -1;
  int bcol = -1;
  //  int snpcol = -1;
  int chrcol = -1;

  vector<string> tokens = tokenizeLine(IN);

  for (int i = 0 ; i < tokens.size() ; i++)
    {
      if ( tokens[i] == "P" )
	pcol = i;
      if ( tokens[i] == "OR" || tokens[i] == "BETA" )
	bcol = i;
//       if ( tokens[i] == "SNP" )
// 	snpcol = i;
      if ( tokens[i] == "CHR" )
	chrcol = i;
      if ( tokens[i] == "BP" )
	bpcol = i;
    }

  int ncol = tokens.size();

  if ( pcol == -1 )
    error("Could not find P field in header");
//   if ( snpcol == -1 )
//     error("Could not find SNP field in header");
  if ( bpcol == -1 )
    error("Could not find BP field in header");
  if ( chrcol == -1 )
    error("Could not find CHR field in header");

  bool no_beta = false;
  if ( bcol == -1 )
    {
      no_beta = true;
      printLOG("Couldn't find OR/BETA field, so reporting all regions\n");
    }

  //  map<int2, string> snpmap;
  map<int2, double> pmap;
  map<int2, double> bmap;

  while ( !IN.eof() )
    {

      vector<string> tokens = tokenizeLine(IN);

      if ( tokens.size() == 0 )
	continue;

      if ( tokens.size() != ncol )
	error("Wrong number of columns in input file");

      double p, b=1;
      int chr, bp;
      //      string snp;

      //snp = tokens[snpcol];
      chr = getChromosomeCode( tokens[chrcol ] );

      if ( ! from_string<double>( p , tokens[pcol] , std::dec ) )
	p = 1;

      if ( ! no_beta )
	{
	  if ( ! from_string<double>( b , tokens[bcol] , std::dec ) )
	    b = 1;
	}

      if ( ! from_string<int>( bp , tokens[bpcol] , std::dec ) )
	error("Problem converting BP value: " + tokens[bpcol] );

      int2 t( chr , bp );

      //      snpmap.insert( make_pair( t , snp ) );
      pmap.insert( make_pair( t , p ) );
      bmap.insert( make_pair( t , b ) );

    }

  IN.close();


  ofstream SUM;
  printLOG("Writing range summary to [ " + par::output_file_name + ".elf.summary ]\n");
  SUM.open( ( par::output_file_name + ".elf.summary").c_str() , ios::out );

  SUM << setw(4) << "CHR" << " "
      << setw(12) << "BP1" << " "
      << setw(12) << "BP2" << " "
      << setw(12) << "BESTP" << "   "
      << "GENES" << "\n";


  map<int2,double>::iterator i = pmap.begin();
  int srange_cnt = 0;
  bool inRange = false;
  Range srange;
  int l = 0;
  int ntot = pmap.size() - 1;
  double bestp = 1;


  while ( i != pmap.end() )
    {

      double pvalue = i->second;

      double coef;
      if ( no_beta )
	coef = 99 ;
      else
	coef = bmap.find( i->first )->second;

      // Look for control enrichment? If so, just flip coef here
      if ( par::rare_test_summary_controls )
	coef = 1 / coef;

      if ( ( ! inRange )
	   && pvalue <= par::rare_test_score_range_threshold
	   && coef > 1 )
	{
	  inRange = true;
	  bestp = pvalue;
	  srange.chr = i->first.p1;
	  srange.start = srange.stop = i->first.p2;
	}
      else if ( inRange
		&& ( pvalue > par::rare_test_score_range_threshold  ||
		     coef < 1 ||
		     i->first.p1 != srange.chr ||
		     l == ntot ) )
	{
	  SUM << setw(4) << srange.chr << " "
	      << setw(12) << srange.start << " "
	      << setw(12) << srange.stop << " "
	      << setw(12) << bestp << "   ";

	  ++srange_cnt;

	  SUM.flush();

	  // Lookup genes in this region?

	  if ( true )
	    {

	      Range r1(srange.chr, srange.start , srange.stop , "dummy");
	      set<Range*> implicated = rangeIntersect(r1,ranges);
	      set<Range*>::iterator ri = implicated.begin();
	      while ( ri != implicated.end() )
		{
		  SUM << (*ri)->name << ",";
		  ++ri;
		}
	    }

	  SUM << "\n";

	  inRange = false;

	  if ( pvalue <= par::rare_test_score_range_threshold
	       && coef > 1 )
	    {
	      srange.chr =  i->first.p1;
	      srange.start = srange.stop = i->first.p2;
	      inRange = true;
	      bestp = pvalue;
	    }
	}




      if ( inRange )
	{
	  srange.stop = i->first.p2;
	  if ( pvalue < bestp && realnum(pvalue) )
	    bestp = pvalue;

	}

      ++l;
      ++i;
    }

  printLOG("Found " + int2str(srange_cnt) + " distinct regions of contiguous association\n");

  SUM.close();

  shutdown();
}

void Plink::elfBaseline()
{

  ofstream ELF;
  string f = par::output_file_name + ".elf.baseline";

  ELF.open(f.c_str(),ios::out);
  ELF.precision(4);

  printLOG("Writing baseline LF SNP count to [ " + f + " ] \n");

  ELF << setw(par::pp_maxfid) << "FID" << " "
      << setw(par::pp_maxiid) << "IID" << " "
      << setw(4) << "CHR" << " "
      << setw(8) << "CNT" << " "
      << setw(8) << "GENO" << " "
      << setw(8) << "RATE" << "\n";

  for (int i = 0 ; i < n ; i++ )
    {
      Individual * person = sample[i];

      int cnt = 0;
      int gcnt = 0;

      map<int,int> chr_cnt;
      map<int,int> chr_gcnt;

      int chr = -1;
      int * p_cnt;
      int * p_gcnt;

      for (int l = 0; l < nl_all; l++)
	{

	  if ( locus[l]->freq > par::rarer_maf_threshold )
	    continue;

	  if ( locus[l]->chr != chr )
	    {
	      chr = locus[l]->chr;
	      chr_cnt.insert( make_pair( chr , 0 ) );
	      chr_gcnt.insert( make_pair( chr , 0 ) );
	      p_cnt = &(chr_cnt.find(chr)->second);
	      p_gcnt = &(chr_gcnt.find(chr)->second);
	    }

	  bool X=false, haploid=false;
	  if ( par::chr_sex[locus[l]->chr] ) X=true;
	  else if ( par::chr_haploid[locus[l]->chr] ) haploid=true;

	  bool s1 = par::SNP_major ? SNP[l]->one[i] : person->one[l];
	  bool s2 = par::SNP_major ? SNP[l]->two[i] : person->two[l];

	  if ( s1 && !s2 )
	    continue;

	  if ( haploid || ( X && person->sex ) )
	    {
	      ++gcnt;
	      ++(*p_gcnt);

	      if ( !s1 )
		{
		  ++cnt;
		  ++(*p_cnt);
		}
	    }
	  else
	    {
	      gcnt +=2;
	      ++(*p_gcnt);
	      ++(*p_gcnt);

	      if ( !s1 )
		{
		  ++cnt;
		  ++(*p_cnt);
		}

	      if ( !s2 )
		{
		  ++cnt;
		  ++(*p_cnt);
		}

	    }
	}

      ELF << setw(par::pp_maxfid) << person->fid << " "
	  << setw(par::pp_maxiid) << person->iid << " "
	  << setw(4) << "G" << " ";
      ELF << setw(8) << cnt << " "
	  << setw(8) << gcnt << " "
	  << setw(8) << (double)cnt / (double)gcnt << "\n";

      map<int,int>::iterator it = chr_cnt.begin();
      while ( it != chr_cnt.end() )
	{
	  int c = it->first;
	  int x = chr_cnt.find( c )->second;
	  int y = chr_gcnt.find( c )->second;

	  ELF << setw(par::pp_maxfid) << person->fid << " "
	      << setw(par::pp_maxiid) << person->iid << " "
	      << setw(4) << c << " ";
	  ELF << setw(8) << x << " "
	      << setw(8) << y << " "
	      << setw(8) << (double)x / (double)y << "\n";

	  ++it;
	}

    }
  ELF.close();

}

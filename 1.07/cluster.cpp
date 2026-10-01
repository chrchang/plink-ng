

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
#include <iomanip>
#include <fstream>
#include <sstream>
#include <set>
#include <algorithm>
#include <cmath>

#include "plink.h"
#include "helper.h"
#include "options.h"
#include "perm.h"
#include "stats.h"

using namespace std;

extern ofstream LOG;

// Helper function: find the maximum distance between two clusters
double cldist(vector<vector<double> > &, vector<int> &, vector<int> &);

// Helper function: group average link
double groupAvgLink(vector<vector<double> > &, vector<int> &, vector<int> &);


// Helper function: are two clusters phenotypically homogeneous?
bool homogeneous_clusters(Plink &, vector<int> &, vector<int> &);

// Do two clusters conform to any --mcc Ncase Ncontrol specification?
bool spec_clusters(Plink &, vector<int> &, vector<int> &);

// Any members of the clusters that can't be matched?
bool pairable_cluster(vector<vector<bool> > & , vector<int>&, vector<int>&);

// Have we already picked somebody from this category?
bool selcon_inds(Plink&, vector<int>&, vector<int>&, set<int>&);


class Neighbour
{
public:
  double dist;
  Individual * neighbour;
  bool operator< (const Neighbour & s2) const
  {  return (dist < s2.dist); }  
};


// Complete-linkage clustering based on average IBS distance

// Extra constraints:
//    --pmerge P     do not merge clusters containing two individuals who differ at this level
//    --mc N         do not let clusters contain more than N individuals
//    --cc  	     do not merge phenotypically identical clusters
//    --mcc N1 N2    do not let cluster contain more than N1 cases and N2 controls
//    --match        external categorical matching criteria
//    --match-type   positive or negative matches
//    --qmatch       external quantitative threshold based
//    --qt           define thresholds for QT matching
//    --pick1        only select one individual from each covariate group
//    --ibm X        identity-by-missingness threshold


void Plink::buildCluster()
{

  ///////////////////////////////////////
  // This is an individual-mode analysis

  if (par::SNP_major) 
    SNP2Ind();


  //////////////////////////////////////
  // Force an initial cluster solution?
   
  // Initially, # of clusters = # of people, unless
  // we are forcing a starting solution
  
  int ni = n;
  
  vector<vector<int> > cl;

  if ( par::force_initial_cluster )
    {
      if (!readClusterFile())
	error("Problem reading --within {file}");

      printLOG("Forcing an initial starting solution from [ " + par::include_cluster_filename + " ]\n");

      set<int> added;
      for (int k=0;k<nk;k++)
	{
	  vector<int> t;
	  for (int i=0;i<ni;i++)
	    {
	      if ( sample[i]->sol == k )
		{
		  t.push_back(i);
		  added.insert(i);
		}
	    }	
	  cl.push_back(t);
	}      
      
      // And now add any remain individuals, in their own clusters,
      // starting from cluster nk onwards

      for (int i=0;i<ni;i++)
	{
	  if ( added.find(i) == added.end() )
	    {
	      vector<int> t(1);
	      t[0] = i;
	      cl.push_back(t);
	    }
	}
      
    }
  else
    {
      for (int i=0;i<ni;i++)
	{
	  vector<int> t(1);
	  t[0] = i;
	  cl.push_back(t);
	}
    }

    
  // T/F matrix (lower diagonal) for whether 
  // two people can be matched, based on p-value
  // constraint and any external criteria
  // pairable[i][j] (default = T)

  vector<vector<bool> > pairable(n);
  for (int i=0; i<n; i++) 
    {
      vector<bool> tmp(n,true);
      pairable[i] = tmp;
    }
  
   
  ///////////////////////////////
  // External matching criteria

  // Determine, in advance, potential pairwise matching

  if (par::bmatch)
  {

    printLOG("Applying categorical matching criteria...\n");

    // Read in each covariate one at a time, 
    // and determine matching
    // we can use the covariate file, as the cluster
    // routine exits after clustering (i.e. so covariates
    // never used)


    // Has the user specified a match-type file? If not, 
    // assume all are positive matches.

    vector<bool> btype(0); 
    if (par::bmatch_usertype)
      {
	checkFileExists(par::bmatch_direction_filename);
	ifstream BT(par::bmatch_direction_filename.c_str(), ios::in);
	while (!BT.eof()) 
	  {
	    string tmp;
	    BT >> tmp;
	    if(BT.eof()) break;	 
	    if (tmp=="+" || tmp=="1") 
	      btype.push_back(true);
	    else
	      btype.push_back(false);
	  }
    BT.close();    
    
    printLOG(int2str(btype.size())+" match-type definitions read from [ "+
	     par::bmatch_direction_filename+" ]\n");
    
      }
    

    int c=0;

    // Swap b-match filename as cluster/within filename
    par::include_cluster_filename = par::bmatch_filename;
    
    while (1)
      { 
	par::mult_clst = ++c;
	if (!readClusterFile()) break;
	
	if (!par::bmatch_usertype)
	  btype.push_back(true);
	
	for (int i=0; i<n-1; i++) 
	  for (int j=i+1; j<n; j++) 
	    {
	      
	      // ->missing means missing on covariate in this context
	      
	      // Simple matching (no usertypes or +-match)
	      if ( btype[c-1] )
		{
		  // +/match
		  if (sample[i]->sol != sample[j]->sol && 
		      (!sample[i]->missing ) && 
		      (!sample[j]->missing ) ) 
		    pairable[i][j] = pairable[j][i] = false;
		}
	      else
		{
		  // -/match
		  if (sample[i]->sol == sample[j]->sol && 
		      (!sample[i]->missing ) && 
		      (!sample[j]->missing ) ) 
		    pairable[i][j] = pairable[j][i] = false;		
		}
	      
	    }
	 
      }
    printLOG("Matched on "+int2str(c-1)+
	     " variables from [ "+par::bmatch_filename+" ]\n");   
  }
 

  if (par::qmatch)
  {

    printLOG("Applying quantitative matching criteria...\n");
	
    vector<double> qt; // number of thresholds specified
     
    checkFileExists(par::qmatch_threshold_filename);
    ifstream QT(par::qmatch_threshold_filename.c_str(), ios::in);
    while (!QT.eof()) 
      {
	double tmp;
	QT >> tmp;
	if(QT.eof()) break;	 
	qt.push_back(tmp);
      }
    QT.close();    
    
    printLOG(int2str(qt.size())+" q-match thresholds read from [ "+
	     par::qmatch_threshold_filename+" ]\n");
    
    
    // Swap q-match filename as covariate file
    par::covar_filename = par::qmatch_filename;
    
    int c=0;   // counter for number of fields in qmatch file

    for (int z=1; z<=qt.size(); z++)
      { 
	par::mult_covar = z;
	if (!readCovariateFile()) break;
	c++;
	
	for (int i=0; i<n-1; i++) 
	  for (int j=i+1; j<n; j++) 
	    {
	      if ( abs( sample[i]->covar - sample[j]->covar ) > qt[c-1] && 
		   (!sample[i]->missing) && 
		   (!sample[j]->missing) ) 
		pairable[i][j] = pairable[j][i] = false;
	    }
	
      }
    
    printLOG("Matched on "+
	     int2str(c)+" quantitative covariates from [ "
	     +par::qmatch_filename +" ]\n");
  }

  if (par::cluster_missing) 
    {
      printLOG("Clustering individuals based on genome-wide IBM\n");
    }
  else
    {
      printLOG("Clustering individuals based on genome-wide IBS\n");
      stringstream s2;
      s2 << "Merge distance p-value constraint = " << par::merge_p << "\n";
      printLOG(s2.str());
    }

  if (par::outlier_detection)
    printLOG("Outlier detection based on neighbours "+int2str(par::min_neighbour)+
	     " to "+int2str(par::max_neighbour)+"\n");  
  
  
  /////////////////////////////////////////////////////////
  // Also, if --pick1 is in effect, we need to read a list from which
  // we can pick only 1 individual
  
  if (par::cluster_selcon)
    {
      // Swap pick1 filename as covariate file
      par::include_cluster_filename = par::cluster_selcon_file;
      par::mult_clst = 1;
      if (!readClusterFile())
	error("Problem reading for --pick1 option");
    }

  // Keep track of what has been selected already
  set<int> selcon;


  /////////////////////////////
  // Set up distance matrices

  // Lower diagonal structure, requires that i > j 
  mdist.resize(n);
  for (int j=0;j<n;j++)
    mdist[j].resize(j);

     
  //////////////////////////////////////////
  // Genome-wide IBS for each pair
  // Either calculate, or re-read from file
  
  vector<double> prop_sig_diff(n);

  // Calculate...
  if (!par::ibd_read)
    {
      int c=0;
      int c2=0;
      
      for (int i1=0; i1<n-1; i1++)
	for (int i2=i1+1; i2<n; i2++)
	  {
	    
	    // Only update message every 100 iterations
	    if (c==c2 || c==np)
	      {
		if (par::cluster_missing)
		  {
		    cout << "IBM calculation: " 
			 << c++ << " of " << np 
			 << "         \r";
		    cout.flush();
		  }
		else
		  {
		    cout << "IBS(g) calculation: " 
			 << c++ << " of " << np 
			 << "         \r";
		    cout.flush();
		  }
		c2+=100;
	      }
	    else
	      ++c;
	    
	    Z IBSg;
	    
	    if (par::cluster_missing) 
	      calcGenomeIBM(sample[i1],sample[i2]);
	    else 
	      {
		// Also calculate IBM as a constraint?
		if (par::cluster_ibm_constraint)
		  {
		    calcGenomeIBM(sample[i1],sample[i2]);
		    if ( dst < par::cluster_ibm_constraint_value )
		      pairable[i1][i2] = pairable[i2][i1] = false;
		  }
		
		// IBS distance (stored in dst)
		IBSg = calcGenomeIBS(sample[i1],sample[i2]);
		
	      }

	    mdist[i2][i1]=dst;

	    //////////////////////////
	    // Is this pair pairable?
	    
	    if (pv < par::merge_p && realnum(pv)) 
	      {   
		// record pair as unpairable
		pairable[i1][i2] = pairable[i2][i1] = false;
		
		// record for both individuals a IBS-based mismatch
		prop_sig_diff[i1]++;
		prop_sig_diff[i2]++;                                
	      }
	  }
      
    }
  else // ... read IBS information from .genome file
    {
      
      

      checkFileExists(par::ibd_file);

      if ( par::ibd_read_minimal )
	printLOG("Reading IBS estimates (minimal format) from [ "
		 +par::ibd_file+" ] \n");
      else
	printLOG("Reading genome-wide IBS estimates from [ "
		 +par::ibd_file+" ] \n");
      

      if ( compressed( par::ibd_file ) )
	par::compress_genome = true;
	   
      ZInput ZINC( par::ibd_file  , par::compress_genome );

      
      map<string,int> mperson;
      for (int i=0; i<n; i++)
	mperson.insert(make_pair( sample[i]->fid+"_"+sample[i]->iid , i )); 
      
      map<Individual*,int> mcode;
      for (int i=0; i<n; i++)
	mcode.insert(make_pair( sample[i] , i )); 

      vector<Individual*> peeps;


      if ( par::ibd_read_minimal )
	{
	  // read in list of people here
	  while ( 1 ) 
	    {
	      vector<string> ids = ZINC.tokenizeLine();
	      if ( ids.size() != 2 ) 
		{
		  string emsg = "Problem with line in [ " + par::ibd_file + " ]\n";
		  for (int i=0;i<ids.size();i++)
		    emsg += ids[i] + " ";
		  error(emsg);
		}

	      string fid = ids[0];
	      string iid = ids[1];

	      if ( fid == "__END" ) 
		break;
	      
	      // Find this person
	      string pcode = fid+"_"+iid;
	      map<string,int>::iterator p = mperson.find(pcode);

	      // Add NULL if this person actually not in 
	      // the current file -- in this case, they 
	      // will be ignored -- but remember we have 
	      // to check for NULLs below and skip those
	      // numbers in that case...

	      if ( p == mperson.end() )
		peeps.push_back( NULL );
	      else
		peeps.push_back( sample[p->second] );
	      
	      // Just in case we have a malformed file
	      
 	      if ( ZINC.endOfFile() )
 		error("Problem with premature stop in file [ " + par::ibd_file + " ]\n");
	      
	    }


	  //////////////////////////////////////////////////////
	  // Now read the actual IBS/PPC values for these peeps

	  if ( peeps.size() != sample.size() )
	    printLOG("Warning -- a different number of people in .genome.min that dataset\n");

	  int size = peeps.size();

	  int p1 = 0, p2 = 1;
	  
	  while ( 1 ) 
	    {
	      double mydst, pv, ibd;

	      vector<string> val = ZINC.tokenizeLine();

	      if ( ZINC.endOfFile() )
		{
		  // Check that p1,p2 counts are as should be...  
		  break;
		}

	      if ( val.size() != 3 ) 
		{
		  string emsg = "Problem with line in [ " + par::ibd_file + " ]\n";
		  for (int i=0;i<val.size();i++)
		    emsg += val[i] + " ";
		  error(emsg);
		}
	      
	      if ( !from_string<double>( mydst, val[0], std::dec ) )
		mydst = 0;

	      if ( !from_string<double>( pv, val[1], std::dec ) )
		pv = 0;

	      if ( !from_string<double>( ibd, val[2], std::dec ) )
		ibd = 0;

	      Individual * person1 = peeps[p1];
	      Individual * person2 = peeps[p2];
	      
	      int pn1 = mcode.find( person1 )->second;
	      int pn2 = mcode.find( person2 )->second;
	      
	      if ( person1 == NULL || person2 == NULL || person1 == person2 ) 
		{
		  // Advance to next pair
		  ++p2;
		  if ( p2 == n )
		    {
		      ++p1;
		      p2=p1+1;
		    }
		  if ( p1==n )
		    break;
		  continue;
		}
	      
// 	      cout << "found " << pn1 << " and " << pn2 << " is " 
// 		   << person1->fid << " " << person1->iid << " x " 
// 		   << person2->fid << " " << person2->iid << "\t"
// 		   << " with " 
// 		   << mydst << " " << pv << "\n";
		
	      // Record IBS distance
	      if ( pn1 > pn2 )
		mdist[pn1][pn2] = mydst;
	      else
		mdist[pn2][pn1] = mydst;
	      
	      
	      //////////////////////////
	      // Is this pair pairable?
	      
	      if (pv < par::merge_p && realnum(pv))
		{   
		  // record pair as unpairable
		  pairable[pn1][pn2] = false;
		  pairable[pn2][pn1] = false;
		  
		  // record for both individuals a IBS-based mismatch
		  prop_sig_diff[pn1]++;
		  prop_sig_diff[pn2]++;                                
		}
	      
	      // Also calculate IBM as a constraint?
	      if (par::cluster_ibm_constraint)
		{
		  
		  calcGenomeIBM(person1,person2);
		  
		  if ( dst < par::cluster_ibm_constraint_value )
		    {
		      pairable[pn1][pn2] = false;
		      pairable[pn2][pn1] = false;
		    }
		}
	      
	      // Advance to next peep-pair
	      ++p2;
	      if ( p2 == n )
		{
		  ++p1;
		  p2=p1+1;
		}

	      // Finished?
	      if ( p1==n )
		break;
	      

	    }
	}
      else
	{
	  
	  // Read in .genome file in verbose mode
	  
	  // We only want FID1,IID1,FID2,IID2 (always first four)
	  // DST and PPC 

	  // Get field codes from header
	  
	  int ppc_code = -1;
	  int dst_code = -1;
	  int col_length = 0;
	  double mydst;

	  vector<string> tokens = ZINC.tokenizeLine();
	  col_length = tokens.size();

	  if ( tokens.size() < 4 || 
	       tokens[0] != "FID1" || 
	       tokens[1] != "IID1" || 
	       tokens[2] != "FID2" || 
	       tokens[3] != "IID2" )
	    error("Problem with header row of .genome file");

	  
	  for ( int i = 4; i<tokens.size(); i++)
	    {
	      if ( tokens[i] == "PPC" )
		ppc_code = i;
	      if ( tokens[i] == "P" )
		ppc_code = i;
	      if ( tokens[i] == "DST" )
		dst_code = i;
	    }
	  
	  if ( ppc_code == -1 || dst_code == -1 )
	    error("Could not find PPC or DST fields in .genome file");
	  
	  // Read each pair at a time
	  while ( ! ZINC.endOfFile() ) 
	    {
	      
	      vector<string> tokens = ZINC.tokenizeLine();

	      if ( tokens.size() == 0 ) 
		continue;

	      if ( col_length != tokens.size() )
		{
		  string strmsg = "";
		  for (int i=0;i<tokens.size();i++)
		    strmsg += tokens[i] + " ";
		  error("Problem reading line in .genome file:\n"+strmsg+"\n");
		}
	      
	      string fid1 = tokens[0];
	      string iid1 = tokens[1];
	      string fid2 = tokens[2];
	      string iid2 = tokens[3];
	      string ipv = tokens[ppc_code];
	      string idst = tokens[dst_code];
	      
	      // Skip any blank rows, or additional header rows
	      if (fid1=="") continue;
	      if (fid1=="FID1") continue;
	      
// 	      if ( ! ( from_string<double>( ibs0 , i0 , std::dec) && 
// 		       from_string<double>( ibs1 , i1 , std::dec) && 
// 		       from_string<double>( ibs2 , i2 , std::dec) ) )	     
// 		{
// 		  error("Problem with line in .genome file, IBS estimates: \n"
// 			+i0+" "+i1+" "+i2+" "+ipv+"\n");
// 		}
	      
	      if ( ! from_string<double>( mydst , idst , std::dec) )
		mydst = 0;
	      
	      if ( ! from_string<double>( pv , ipv , std::dec) )
		pv = 1;
	      
	      // Calculate proportion IBS matching

// 	      if (par::cluster_euclidean)
// 		mydst = (ibs2*2+ibs1*0.5)/(ibs2*2+ibs1+ibs0);
// 	      else 
// 		mydst = (ibs2+ibs1*0.5)/(ibs2+ibs1+ibs0);
	      
	      
	      map<string,int>::iterator person1 = mperson.find(fid1+"_"+iid1);
	      map<string,int>::iterator person2 = mperson.find(fid2+"_"+iid2);
	      
	      if ( person1 == mperson.end() || person2 == mperson.end() || person1 == person2 ) 
		continue;
	      
	      // Record IBS distance
	      if ( person1->second > person2->second )
		mdist[person1->second][person2->second] = mydst;
	      else
		mdist[person2->second][person1->second] = mydst;
	      
	      
	      //////////////////////////
	      // Is this pair pairable?
	      
	      if (pv < par::merge_p && pv==pv) 
		{   
		  // record pair as unpairable
		  pairable[person1->second][person2->second] = false;
		  pairable[person2->second][person1->second] = false;
		  
		  // record for both individuals a IBS-based mismatch
		  prop_sig_diff[person1->second]++;
		  prop_sig_diff[person2->second]++;                                
		}
	      
	      // Also calculate IBM as a constraint?
	      if (par::cluster_ibm_constraint)
		{
		  
		  calcGenomeIBM(sample[person1->second],sample[person2->second]);
		  
		  if ( dst < par::cluster_ibm_constraint_value )
		    {
		      pairable[person1->second][person2->second] = false; 
		      pairable[person2->second][person1->second] = false;
		    }
		}
	      
	    } // Read next line in .genome
	}
      
      
      ZINC.close();
      
      
      /////////////////////////////////////////////
      // Check that every pair in the dataset has 
      // actually been assigned a value -- i.e. check 
      // for 0 IBS codes, etc.

      
    }
  

  ///////////////////////////////////
  // IBS permutation test

  if ( par::ibs_test ) 
    {
      // If we were called by permutationIBSTest(), 
      // now it is time to return 
      
      return;
    }


  ///////////////////////////////////
  // Display matrix of IBS distances
  
  if (par::matrix)
    {
      string f;
      
      if (par::cluster_missing) 
	f = par::output_file_name+ ".mdist.missing";
      else if (par::distance_matrix)
	f = par::output_file_name+ ".mdist";
      else
	f = par::output_file_name+ ".mibs";

      if (!par::cluster_missing)
	{
	  if (par::distance_matrix)
	    printLOG("Writing IBS distance matrix to [ "+f + " ]\n");
	  else
	    printLOG("Writing IBS similarity matrix to [ "+f + " ]\n");
	}
      else
	printLOG("Writing IBM distance matrix to [ "+f + " ]\n");

      ofstream MAT(f.c_str(),ios::out);
      MAT.clear();
      
      for (int i=0;i<mdist.size();i++)
	{
	  for (int j=0;j<mdist.size();j++)
	    {
	      if ( par::distance_matrix ) 
		{
		  // Distances
		  if (i>j)
		    MAT << 1 - mdist[i][j] << " ";
		  else if (i==j)
		    MAT << 0 << " ";
		  else 
		    MAT << 1 - mdist[j][i] << " ";
		}
	      else
		{
		  // Similarities
		  if (i>j)
		    MAT << mdist[i][j] << " ";
		  else if (i==j)
		    MAT << 1 << " ";
		  else 
		    MAT << mdist[j][i] << " ";
		}
	    }
	  MAT << "\n";
	}
      MAT.close();

    }


  ////////////////////////////////////
  // Determine how many pairable pairs 
  // we have now
  
  if (!par::cluster_missing)
    {
      int paircount = 0;
      for (int i=0; i<n-1; i++) 
	for (int j=i+1; j<n; j++) 
	  if (pairable[i][j]) paircount++;
      printLOG("Of these, "+int2str(paircount)+" are pairable based on constraints\n");
    }
 


  //////////////////////////
  // Outlier detection
  
  if (par::outlier_detection)
    {
      printLOG("Writing individual neighbour/outlier statatistics to [ " +
	       par::output_file_name + ".nearest ]\n");
      
      vector<vector<double> > min_dst(n);
      vector<vector<double> > zmin_dst(n);
      vector<vector<Individual*> > min_ind(n);
      
      if (par::max_neighbour > n-1)
	error("Nearest neighbour range specified as [ "+int2str(par::max_neighbour)
	      +" ] but only [ "+int2str(n)+" ] individuals in sample.");
      
      for (int k=par::min_neighbour;k<=par::max_neighbour;k++)
	{  
	  
	  // Consider each person
	  for (int i=0;i<n;i++)
	    {
	      
	      vector<Neighbour> ibs(n-1);
	      
	      int c=0;
	      for (int j=0;j<n;j++)
		if (i!=j) 
		  {
		    if ( i>j )
		      ibs[c].dist = mdist[i][j];
		    else
		      ibs[c].dist = mdist[j][i];
		    
		    ibs[c].neighbour = sample[j];
		    c++;
		  }

	      sort(ibs.begin(),ibs.end());
	      min_dst[i].push_back(ibs[ibs.size() - k].dist);
	      min_ind[i].push_back(ibs[ibs.size() - k].neighbour);
	    }
	  
	  // Calculate mean and variance of min_dst to 
	  // give Z-scores
	  
	  double mean = 0;
	  double var = 0;
	  for (int i=0; i<n; i++)
	    mean += min_dst[i][min_dst[i].size()-1];
	  mean /= (double)n;
	  for (int i=0; i<n; i++)
	    var += (min_dst[i][min_dst[i].size()-1]-mean)*(min_dst[i][min_dst[i].size()-1]-mean);
	  var /= (double)(n-1);
	  for (int i=0; i<n; i++)
	    zmin_dst[i].push_back( ( min_dst[i][min_dst[i].size()-1] - mean ) / sqrt(var) ) ;
	  
	}
      
      // Second measure: based on significance test
      
      // Proportion of rest of sample with whom significant difference at 'pmerge' threshold
      // pv might be NaN, but only if very small # of markers is used -- ignore, as
      // all values will be meaningless in any case
      
      if (!par::cluster_missing)
	for (int i=0;i<n;i++)
	  prop_sig_diff[i] /= (double)(n-1);
	
      
      // And output to a file
      
      ofstream MD((par::output_file_name+".nearest").c_str(),ios::out);
      MD.clear();
      MD.precision(4);
      
      MD << setw(12) << "FID" << " "
	 << setw(12) << "IID" << " "
	 << setw(6)  << "NN" << " "
	 << setw(12) << "MIN_DST" << " "
	 << setw(12) << "Z"  << " "
	 << setw(12) << "FID2" << " "
	 << setw(12) << "IID2" << " ";
      if (!par::cluster_missing)
	MD << setw(12) << "PROP_DIFF" << " ";
      MD << "\n";
      
      for (int i=0; i<n; i++)
	for (int k=0;k<min_dst[0].size();k++)
	  {
	    MD << setw(12) << sample[i]->fid  << " "
	       << setw(12) << sample[i]->iid  << " "
	       << setw(6) << par::min_neighbour+k << " "
	       << setw(12) << min_dst[i][k]  << " "
	       << setw(12) << zmin_dst[i][k]  << " "
	       << setw(12) << min_ind[i][k]->fid << " "
	       << setw(12) << min_ind[i][k]->iid << " ";
	    if (!par::cluster_missing)
	      MD << setw(12) << prop_sig_diff[i] << " ";
	    MD << "\n";
	  }
      MD.close();
    }



  //////////////////////////
  // Cluster analysis
  
  if ( par::cluster ) 
    {

      int c=1;

      bool done=false;
      
      // Matrix of solutions
      vector< vector<int> > sol(ni);
      for (int i=0;i<ni;i++) sol[i].resize(ni);
      
      vector<double> hist(1);
      
      // Build solution
      for (int i=0; i<cl.size(); i++)
	for (int j=0; j<cl[i].size(); j++)
	  sol[cl[i][j]][0] = i;
      
      
      printLOG("Writing cluster progress to [ "+par::output_file_name + ".cluster0 ]\n");
      ofstream CLST((par::output_file_name+".cluster0").c_str(),ios::out);
      CLST.clear();
      
      
      while(!done)
	{
	  
	  double dmin = -999;
	  
	  int imin=-1;
	  int jmin=-1;
	  
	  
	  // 1. Find min/max distance between pairable clusters
	  
	  for (int i=0; i<cl.size()-1; i++)
	    for (int j=i+1; j<cl.size(); j++)
	      {	
		
		// Cluster on IBS: group average link or complete linkage?
		double d = par::cluster_group_avg ? groupAvgLink(mdist,cl[i],cl[j]) : cldist(mdist,cl[i],cl[j]);
		
		// Are these individuals/clusters more similar AND pairable?
		
		if ( d>dmin && pairable_cluster(pairable,cl[i],cl[j])  )
		  {
		    
		    // And will the max cluster size requirement be fulfilled?
		    
		    if (par::max_cluster_size==0 ||  
			(( cl[i].size()+cl[j].size()) <= par::max_cluster_size) )
		      {
			
			// And will the basic phenotypic matching requirement be fulfilled?
			
			if ( (!par::cluster_on_phenotype) 
			     || (!homogeneous_clusters((*this),cl[i],cl[j]))) 
			  {
			    
			    // What about the --mcc clustering 
			    if ( (!par::cluster_on_mcc) || spec_clusters( (*this),cl[i],cl[j]) )
			      {
				
				// And what about pick1 constrains? (this must be final constraint)
				if ( (!par::cluster_selcon) || selcon_inds( (*this),cl[i],cl[j],selcon))
				  {
				    imin=i;
				    jmin=j;
				    dmin=d;
				  }
			      }
			  }
		      }
		    
		  }
		
	      }
	  
	  // Did we get a merge?
	  if (imin==-1) {
	    done=true;
	    //printLOG("Cannot make clusters that satisfy constraints at step "+int2str(c)+"\n");	
	    goto done_making_clusters;
	  }
	  
	  // Save merge distance 
	  hist.push_back(dmin);
	  
	  // Add to list of selected categories
	  if (par::cluster_selcon)
	    {
	      if (cl[imin].size() == 1 )
		selcon.insert( sample[cl[imin][0]]->sol );
	      if (cl[jmin].size() == 1 )
		selcon.insert( sample[cl[jmin][0]]->sol );
	    }
	  
	  // 2. Join these clusters 
	  for(int j=0;j<cl[jmin].size();j++)
	    cl[imin].push_back(cl[jmin][j]);
	  cl.erase(cl.begin()+jmin);
	  if (cl.size()==1 || cl.size()==par::max_cluster_N) done=true;
	  
	  
	  
	  // List entire sample
	  CLST << "Merge step " << c << "\t" << hist[c];
	  
	  // Build solution
	  for (int i=0; i<cl.size(); i++)
	    for (int j=0; j<cl[i].size(); j++)
	      {
		sol[cl[i][j]][c] = i;
	      }
	  
	  
	  // Calculate average within/between cluster distances
	  double between = 0, within = 0; 
	  int withinN = 0, betweenN = 0;
	  
	  for (int j1=0; j1<sol.size(); j1++)
	    for (int j2=0; j2<sol.size(); j2++)
	      {
		if (j1 < j2)
		  {
		    if(sol[j1][c] == sol[j2][c])
		      {
			within += mdist[j2][j1];
			withinN++;		  
		      }
		    else
		      {
			between += mdist[j2][j1];
			betweenN++;		  
		      }
		  }
	      }
	  
	  CLST << "\t" << between/(double)betweenN 
	       << "\t" << within/(double)withinN 
	       << "\t" << ( between/(double)betweenN )  / ( within/(double)withinN )
	       << "\n";
	  
      
	  // Next merge
	  c++;
	}
      
      
    done_making_clusters:
      
      
      CLST.close();
      
      
      //////////////////////////////////
      // Best solution is final solution
      
      int best = hist.size()-1;
      
      if (!par::cluster_missing)
	{
	  printLOG("Writing cluster solution (1) [ " 
		   + par::output_file_name + ".cluster1 ]\n");
	  CLST.open((par::output_file_name+".cluster1").c_str(),ios::out);
	  CLST.clear();
	  
	  for (int i=0; i<cl.size(); i++)
	    {
	      CLST << "SOL-" << i << "\t"; 
	      for (int j=0; j<cl[i].size(); j++)
		{
		  CLST << " " 
		       << sample[cl[i][j]]->fid << "_"
		       << sample[cl[i][j]]->iid;
		  if (par::cluster_on_phenotype || par::cluster_on_mcc)
		    CLST << "("
			 << (int)sample[cl[i][j]]->phenotype
			 << ")";
		  
		}
	      CLST << "\n";
	    }
	  
	  CLST.close();
	  
	  
	  printLOG("Writing cluster solution (2) [ " 
		   + par::output_file_name + ".cluster2 ]\n");
	  CLST.open((par::output_file_name+".cluster2").c_str(),ios::out);
	  
	  CLST.clear();
	  
	  for (int j=0; j<sol.size(); j++)
	    {
	      // Display...
	      CLST << sample[j]->fid << " "
		   << sample[j]->iid << "\t" 
		   << sol[j][best] << "\n";

	      // Keep track of this (might be needed if MDS plot done)
	      sample[j]->sol = sol[j][best];
	    }
	  
	  CLST.close();
	}


  
      if (!par::cluster_missing)
	{
	  printLOG("Writing cluster solution (3) [ " 
		   + par::output_file_name + ".cluster3 ]\n");
	  CLST.open((par::output_file_name+".cluster3").c_str(),ios::out);
	}
      else
	{
	  printLOG("Writing cluster solution (3) [ " 
		   + par::output_file_name + ".cluster3.missing ]\n");
	  CLST.open((par::output_file_name+".cluster3.missing").c_str(),ios::out);
	}
      CLST.clear();
      
      for (int j=0; j<sol.size(); j++)
	{
	  // Display...
	  CLST << sample[j]->fid << " "
	       << sample[j]->iid << "\t";
	  
	  for (int i=0; i<sol[0].size(); i++)
	    CLST << sol[j][i] << " ";
	  
	  CLST << "\n";
	}
      CLST << "\n";
      
      CLST.close();
      
      
    }
  

  //////////////////////////////////////////////////////////
  // Produce MDS plot?
  
  if ( par::cluster_plot )
    generateMDS();
    
  
  // Shutdown now
  shutdown();
  
}


double cldist(vector<vector<double> > & d, 
	      vector<int> & a, 
	      vector<int> & b)
{
  // Compare based on first metric, but also return paired second
  double l;
  l = a[0]>b[0] ? d[a[0]][b[0]] : d[b[0]][a[0]]; 
  
  for (int i=0; i<a.size(); i++)
    for (int j=0; j<b.size(); j++)
      {

	if ( a[i] > b[j] )
	  {
	    if ( d[a[i]][b[j]] < l ) l = d[a[i]][b[j]];
	  }
	else
	  {
	    if ( d[b[j]][a[i]] < l ) l = d[b[j]][a[i]];	    
	  }

      }
  return l;
}


double groupAvgLink(vector<vector<double> > & d, 
		    vector<int> & a, 
		    vector<int> & b)
{
  
  double s = 0;
  
  for (int i=0; i<a.size(); i++)
    for (int j=0; j<b.size(); j++)
      {
	
	if ( a[i] > b[j] )
	  {
	    s += d[a[i]][b[j]];
	  }
	else
	  {
	    s += d[b[j]][a[i]];	    
	  }
	
      }
  
  return 1.0 / ( a.size() * b.size() ) *  s ;
}



bool homogeneous_clusters(Plink & P, vector<int> & a, vector<int> & b)
{

  // Determine how to handle missing phenotypes?

  bool homogeneous = true;
  for (int i=0; i<a.size(); i++)
    for (int j=0; j<b.size(); j++)
      {
	if ( (P.sample[a[i]]->phenotype != P.sample[b[j]]->phenotype) 
	     && (!P.sample[a[i]]->missing) && (!P.sample[b[j]]->missing) )
	  homogeneous = false;
      }
  return homogeneous;
}


bool spec_clusters(Plink & P, vector<int> & a, vector<int> & b)
{
  
  // Missing individuals will be treated as unaffected

  int ncase = 0, ncontrol = 0;

  for (int i=0; i<a.size(); i++)
    if (P.sample[a[i]]->aff) ncase++;
    else ncontrol++;
  
  for (int j=0; j<b.size(); j++)
    if (P.sample[b[j]]->aff) ncase++;
    else ncontrol++;
  
  if (ncase <= par::max_cluster_case && 
      ncontrol <= par::max_cluster_control) 
    return true;
  else
    return false;
}


bool pairable_cluster(vector<vector<bool> > & pairable, vector<int> & a, vector<int> & b)
{
  for (int i=0; i<a.size(); i++)
    for (int j=0; j<b.size(); j++)
       if (!pairable[a[i]][b[j]]) return false; 
  return true;
}

bool selcon_inds(Plink & P, vector<int> & a, vector<int> & b, set<int> & inc)
{
  // Only need to check for singletons (i.e. once somebody is in a cluster, 
  // they must have already passed this test)

  if ( a.size() == 1 )
    {
      // Individual already in?
      if ( inc.find(P.sample[a[0]]->sol) != inc.end() )
	return false;
    }
  else if ( b.size() == 1 )
    {
      if ( inc.find(P.sample[b[0]]->sol) != inc.end() )
	return false;
    }
  
  if ( a.size() == 1 && b.size() == 1 ) 
    if ( P.sample[a[0]]->sol == P.sample[b[0]]->sol )
      return false;
  
  return true;
}




void Plink::permutationIBSTest(Perm & perm)
{
  // Take the IBS distance matrix, and ask (by permutation)
  // where the average difference between two groups is larger
  // than we would expect by chance
  
  // i.e. statistic = average between group IBS distance
  //      permutation = label swapping
  //      1-sided test, asking whether people between
  //      groups are *less* similar than we'd expect

  /////////////////////////////////
  // Calculate distances
  // (will exit before clustering)

  buildCluster();


  ////////////////
  // Perform test
  
  perm.setTests(12);
  perm.setPermClusters(*this);
  perm.originalOrder();

  // Tests (1 sided), where ">" means less similar? 
  // as tests are based on 1-f(mdist)

  // 0  1a. Case/control <  all others
  // 1  1b. Case/control >  all others
 
  // 2  2a. Case/case < control/control
  // 3  2b. Case/case > control/control
 
  // 4  3a. Case/case < all others
  // 5  3b. Case/case > all others

  // 6  4a. Control/control < all others
  // 7  4b. Control/control > all others

  // 8  5a. Case/case < Case/control
  // 9  5b. Case/case > Case/control

  // 10  6a. Control/control < Case/control
  // 11  6b. Control/control > Case/control

  ///////////////////////////
  // Original test statistics

  vector<double> original(12,0);
  
  double bg_mean = 0; 
  double ig1_mean = 0;
  double ig2_mean = 0;
  
  int bg_n = 0;
  int ig1_n = 0;
  int ig2_n = 0;

  double bg_var = 0; 
  double ig1_var = 0;
  double ig2_var = 0;


  // Add addition 11 tests here...

  for (int i=0; i<n-1; i++)
      for (int j=i+1; j<n; j++)
      {
	  if ( sample[i]->aff != sample[j]->aff )
	  {
	      bg_mean += mdist[j][i];
	      bg_n++; 
	  }
	  else if ( sample[i]->aff )
	  {
	      ig2_mean += mdist[j][i];
	      ig2_n++; 
	  }
	  else
	  {
	      ig1_mean += mdist[j][i];
	      ig1_n++; 
	  }
      }
  
  // 0  1a. Case/control <  all others
  // 1  1b. Case/control >  all others
 
  // 2  2a. Case/case < control/control
  // 3  2b. Case/case > control/control
 
  // 4  3a. Case/case < all others
  // 5  3b. Case/case > all others

  // 6  4a. Control/control < all others
  // 7  4b. Control/control > all others

  // 8  5a. Case/case < Case/control
  // 9  5b. Case/case > Case/control

  // 10  6a. Control/control < Case/control
  // 11  6b. Control/control > Case/control

  original[0] -= bg_mean;
  original[1] += bg_mean;

  original[2] += ig1_mean - ig2_mean;
  original[3] += ig2_mean - ig1_mean;

  original[4] -= ig2_mean;
  original[5] += ig2_mean;
  
  original[6] -= ig1_mean;
  original[7] += ig1_mean;

  original[8] += bg_mean - ig2_mean;
  original[9] += ig2_mean - bg_mean;

  original[10] += bg_mean - ig1_mean;
  original[11] += ig1_mean - bg_mean;


  if (bg_n==0)
    error("No between group individuals observed");

  double tot_mean = (bg_mean+ig1_mean+ig2_mean)
    /(double)(bg_n+ig1_n+ig2_n);
  
  bg_mean /= (double)bg_n;
  ig1_mean /= (double)ig1_n;
  ig2_mean /= (double)ig2_n;
  
  for (int i=0; i<n-1; i++)
    for (int j=i+1; j<n; j++)
      {
	if ( sample[i]->aff != sample[j]->aff )
	  bg_var += ( mdist[j][i] - bg_mean ) * ( mdist[j][i] - bg_mean );
	else if ( sample[i]->aff )
	  ig2_var += ( mdist[j][i] - ig2_mean ) * ( mdist[j][i] - ig2_mean );
	else
	  ig1_var += ( mdist[j][i] - ig1_mean ) * ( mdist[j][i] - ig1_mean );
      }
  

  // Total sum of squares
  double total_ss = bg_var + ig2_var + ig1_var;
  
  // Between sum of squares
  double ig_mean = (ig1_mean * ig1_n + ig2_mean * ig2_n) 
    / ( double ) ( ig1_n + ig2_n );
  
  double between_ss = (double)bg_n * ( bg_mean - tot_mean ) * ( bg_mean - tot_mean ) +
    (double)(ig1_n+ig2_n) * ( ig_mean - tot_mean ) * ( ig_mean - tot_mean );
  
  bg_var /= (double)(bg_n-1);
  ig2_var /= (double)(ig2_n-1);
  ig1_var /= (double)(ig1_n-1);
  
  printLOG("\nBetween-group IBS (mean, SD) = "
	   +dbl2str(bg_mean)+", "+dbl2str(sqrt(bg_var))+"\n");
  printLOG("In-group (2) IBS (mean, SD) = "
	   +dbl2str(ig2_mean)+", "+dbl2str(sqrt(ig2_var))+"\n");
  printLOG("In-group (1) IBS (mean, SD) = "
	   +dbl2str(ig1_mean)+", "+dbl2str(sqrt(ig1_var))+"\n");
  printLOG("Approximate proportion of variance between group = " 
	   +dbl2str(between_ss / total_ss)+"\n");
  
  
  ////////////////////
  // Begin permutations
  
  bool finished = false;

  while(!finished)
  {
      
      vector<double> pr(12,0);
      
      // Permute
      perm.permuteInCluster();
      
      // Retest
      bg_mean = ig1_mean = ig2_mean = 0;
      
      for (int i=0; i<n-1; i++)
	  for (int j=i+1; j<n; j++)
	  {
	      if ( sample[i]->pperson->aff != sample[j]->pperson->aff )
	      {
		  bg_mean += mdist[j][i];
	      }
	      else if ( sample[i]->pperson->aff )
	      {
		  ig2_mean += mdist[j][i];
	      }
	      else
	      {
		  ig1_mean += mdist[j][i];
	      }
	  }
      
      pr[0] -= bg_mean;  
      pr[1] += bg_mean;  // are case/control more similar?
      
      pr[2] += ig1_mean - ig2_mean;
      pr[3] += ig2_mean - ig1_mean;
      
      pr[4] -= ig2_mean;
      pr[5] += ig2_mean;
      
      pr[6] -= ig1_mean;
      pr[7] += ig1_mean;
      
      pr[8] += bg_mean - ig2_mean;
      pr[9] += ig2_mean - bg_mean;
      
      pr[10] += bg_mean - ig1_mean;    
      pr[11] += ig1_mean - bg_mean;


      ////////////////////////////////
      // Standard permutation counting
      
      finished = perm.update(pr,original);
      
    }
  
  if (!par::silent)
    cout << "\n\n";
  
  
  ////////////////////////////
  // Display permuted p-values

  printLOG("IBS group-difference empirical p-values:\n\n");
  printLOG(" T1: Case/control less similar                     p = " + dbl2str(perm.pvalue(0)) +"\n");
  printLOG(" T2: Case/control more similar                     p = " + dbl2str(perm.pvalue(1)) +"\n\n");

  printLOG(" T3: Case/case less similar than control/control   p = " + dbl2str(perm.pvalue(2)) +"\n");
  printLOG(" T4: Case/case more similar than control/control   p = " + dbl2str(perm.pvalue(3)) +"\n\n");

  printLOG(" T5: Case/case less similar                        p = " + dbl2str(perm.pvalue(4)) +"\n");
  printLOG(" T6: Case/case more similar                        p = " + dbl2str(perm.pvalue(5)) +"\n\n");

  printLOG(" T7: Control/control less similar                  p = " + dbl2str(perm.pvalue(6)) +"\n");
  printLOG(" T8: Control/control more similar                  p = " + dbl2str(perm.pvalue(7)) +"\n\n");

  printLOG(" T9: Case/case less similar than case/control      p = " +dbl2str(perm.pvalue(8)) +"\n" );
  printLOG("T10: Case/case more similar than case/control      p = " +dbl2str(perm.pvalue(9))  +"\n\n");

  printLOG("T11: Control/control less similar than case/control p = " + dbl2str(perm.pvalue(10)) +"\n");
  printLOG("T12: Control/control more similar than case/control p = " + dbl2str(perm.pvalue(11)) +"\n");

}
  

void Plink::groupGenome()
{

  // Read from a (non-verbose) genome file

  checkFileExists(par::ibd_file);
  
  if ( par::ibd_read_minimal )
    printLOG("Reading IBS estimates (minimal format) from [ "
	     +par::ibd_file+" ] \n");
  else
    printLOG("Reading genome-wide IBS estimates from [ "
	     +par::ibd_file+" ] \n");
  
  ifstream INC;
  INC.open(par::ibd_file.c_str());
  

  map<string,int> mperson;
  for (int i=0; i<n; i++)
    mperson.insert(make_pair( sample[i]->fid+"_"+sample[i]->iid , i )); 
      
  map<Individual*,int> mcode;
  for (int i=0; i<n; i++)
    mcode.insert(make_pair( sample[i] , i )); 
  
  vector<Individual*> peeps;
  
  // We wish to read in an NxN matrix, and convert it to a KxK one
  
  matrix_t dk( nk );
  sizeMatrix(dk,nk,0);
  for (int j=0; j<nk; j++)
    dk[j].resize(j,0);

  table_t dkn;
  sizeTable(dkn,nk,0);
  for (int j=0; j<nk; j++)
    dkn[j].resize(j,0);


  // Read in .genome file in verbose mode
	  
  // We only want FID1,IID1,FID2,IID2 (always first four)
  // DST and PPC 

  // Get field codes from header
  
  int dst_code = -1;
  int col_length = 0;
  double mydst;
  
  vector<string> tokens = tokenizeLine(INC);
  col_length = tokens.size();
  
  if ( tokens.size() < 4 || 
       tokens[0] != "FID1" || 
       tokens[1] != "IID1" || 
       tokens[2] != "FID2" || 
       tokens[3] != "IID2" )
    error("Problem with header row of .genome file");

	  
  for ( int i = 4; i<tokens.size(); i++)
    {
      if ( tokens[i] == "DST" )
	dst_code = i;
    }
  
  if ( dst_code == -1 )
    error("Could not find DST fields in .genome file");
  
  // Read each pair at a time
  while ( ! INC.eof() ) 
    {
      
      vector<string> tokens = tokenizeLine(INC);
      
      if ( tokens.size() == 0 ) 
	continue;
      
      if ( col_length != tokens.size() )
	{
	  string strmsg = "";
	  for (int i=0;i<tokens.size();i++)
	    strmsg += tokens[i] + " ";
	  error("Problem reading line in .genome file:\n"+strmsg+"\n");
	}
      
      string fid1 = tokens[0];
      string iid1 = tokens[1];
      string fid2 = tokens[2];
      string iid2 = tokens[3];
      string idst = tokens[dst_code];
      
      if (fid1=="") continue;
      
      if ( ! from_string<double>( mydst , idst , std::dec) )
	mydst = 0;
      
      map<string,int>::iterator person1 = mperson.find(fid1+"_"+iid1);
      map<string,int>::iterator person2 = mperson.find(fid2+"_"+iid2);
      
      if ( person1 == mperson.end() || 
	   person2 == mperson.end() || 
	   person1 == person2 ) 
	continue;
      
      int k1 = sample[ person1->second ]->sol;
      int k2 = sample[ person2->second ]->sol;

      if ( k1 < 0 || k2 < 0 || k1 == k2 ) 
	continue;
      
      if ( k2 > k1 ) 
	{
	  int tmp = k2;
	  k2 = k1;
	  k1 = tmp;
	}
      
      // Record IBS distance
      
       dk[k1][k2] = mydst;
       ++dkn[k1][k2];
      
    } // Read next line in .genome
  
  INC.close();

  for (int i=0; i<nk; i++)
    for (int j=0; j<i; j++)
      {
	if ( dkn[i][j] > 0 )
	  dk[i][j] /= (double)dkn[i][j];
      }
  
  // Output a dummy .genome file

  ofstream GOUT0;
  GOUT0.open( (par::output_file_name + ".plst").c_str(), ios::out);
  printLOG("Writing person include list to [ " + par::output_file_name + ".plst ]\n");
  
  ofstream GOUT1;
  GOUT1.open( (par::output_file_name + ".clst").c_str(), ios::out);
  printLOG("Writing cluster list to [ " + par::output_file_name + ".clst ]\n");


  map<int,int> k2i;
  for (int i=0;i<n; i++)
    {
      int j = sample[i]->sol;
      if ( j < 0 ) 
	continue;
      if ( k2i.find(j) == k2i.end() )
	{
	  k2i.insert(make_pair(j,i));
	  GOUT0 << sample[i]->fid << " " 
		<< sample[i]->iid << "\n";
	  GOUT1 << sample[i]->fid << " " 
		<< sample[i]->iid << " "
		<< kname[j] << "\n";
	}
    }
  GOUT0.close();
  GOUT1.close();

  ofstream GOUT;
  GOUT.open( (par::output_file_name + ".genome").c_str(), ios::out);
  printLOG("Writing grouped .genome file to [ " + par::output_file_name + ".genome ]\n");
  GOUT << setw(par::pp_maxfid) << "FID1" << " "
       << setw(par::pp_maxiid) << "IID1" << " "
       << setw(par::pp_maxfid) << "FID2" << " "
       << setw(par::pp_maxiid) << "IID2" << " "
       << setw(8) << "DST" << " "
       << setw(8) << "PPC" << "\n";
  
  for (int i=0; i<nk; i++)
    for (int j=0; j<i; j++)
      {
	Individual * s1 = sample[k2i.find(i)->second];
	Individual * s2 = sample[k2i.find(j)->second];

	GOUT << s1->fid << " " << s1->iid << " " 
	     << s2->fid << " " << s2->iid << " ";

// 	cout << i << " " << j << " ";
// 	cout << dk.size() << " " << dk[i].size() << "\n";

// 	cout << dkn.size() << " " << dkn[i].size() << "\n";
//	cout << dk[i][j] << " of " << dkn[i][j] << "\n";

	GOUT << dk[i][j] << " 1\n";
      }
  
  GOUT.close();

}

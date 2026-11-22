

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
#include <string>
#include <cstdlib>
#include <cmath>

#include "options.h"
#include "helper.h"
#include "plink.h"
#include "phase.h"
#include "stats.h"
#include "zed.h"
#include "model.h"
#include "logistic.h"
#include "linear.h"

extern Plink * PP;



class Var {
public:
  string snp;
  string a1;
  bool operator< (const Var & b) const
  {
    return (snp < b.snp);
  }
  bool operator== (const Var & b) const
  {
    return (snp == b.snp && a1 == b.a1 );
  }
};


// Helper functions

void setUpQScoring( map<string,double> & ,
		   vector<double2> &,
		   vector<string> &);



void Plink::processDosageFile()
{



  // For scroring procedure, if used
  
  map<Var, vector<int> >  qflag;
  map<Var, double> wt;
  matrix_t scores;

  map<string,double> qscore;
  vector<double2> qthresh;
  vector<string> qlabel;

  if ( par::score_risk )
    {

      if ( par::score_risk_on_qrange )
	setUpQScoring( qscore, qthresh, qlabel);      
	
      
      // Read in score(s)
      ZInput zin( par::score_risk_file , compressed( par::score_risk_file ) );
      
      while ( ! zin.endOfFile() )
	{
	  vector<string> tok = zin.tokenizeLine();
	  if ( tok.size() != 3 )
	    continue;

	  Var v;
	  v.snp = tok[0];
	  v.a1 = tok[1];

	  double s;
	  if ( ! from_string<double>( s, tok[2] , std::dec ) )
	    continue;
	  
	  // Store weight
	  wt.insert(make_pair( v, s ));
	  
	  // Which scores should we place this in?
	  if ( par::score_risk_on_qrange )
	    {
	      vector<int> qq;
	      
	      // Can we find a q-score?
	      map<string,double>::iterator i1 = qscore.find( v.snp ); 
	      if ( i1 == qscore.end() ) 
		continue;
	      double sc = i1->second;
	      
	      for ( int q = 0 ; q < qthresh.size() ; q++)
		{
		  if ( sc >= qthresh[q].p1 && sc <= qthresh[q].p2 )
		    qq.push_back(q);
		}
	      qflag.insert(make_pair(v,qq));       
	    }
	  // Read next score
	}            
	zin.close();
			
	printLOG("Done reading scores\n");
	
    }
  
  
  printLOG("\nReading dosage information from [ " + par::dosage_file + " ]\n");

  // Expect a single file, with header
  
  bool compressed_in = false;
  bool compressed_out = false;
  bool filelist = false;
  bool countOccur = false;
  bool header = true;
  bool sepHeader = false;
  int skip0 = 0;
  int skip1 = 0;
  int skip2 = 0;
  bool dosageScale2 = true;
  bool snpBatch = false;


  // {skip0} SNP {skip1} A1 A2 {skip2} DATA... 
  
  OptionSet * dosage_opt = par::opt.getOptions("DOSAGE");

  if ( dosage_opt->isSet("Z") ) 
    compressed_in = compressed_out = true;
  
  if ( dosage_opt->isSet("Zout") ) 
    {
      compressed_out = true;
    }
  if ( dosage_opt->isSet("Zin") ) 
    compressed_in = true;
  
  if ( dosage_opt->isSet("list") )
    filelist = true;
  
  if ( dosage_opt->isSet("occur") )
    countOccur = true;
  
  if ( dosage_opt->isSet("noheader") )
    header = false;
  
  if ( dosage_opt->isSet("sepheader") )
    sepHeader = true;

  if ( sepHeader && ! filelist )
    error("The 'sepheader' option requires the 'list' option to be set");
  
  if ( sepHeader && ! header )
    error("Cannot specify both 'sepheader' and 'noheader'");

  if ( dosage_opt->getValue("skip0") != "" ) 
    skip0 = atoi( dosage_opt->getValue("skip0").c_str() );
  if ( dosage_opt->getValue("skip1") != "" ) 
    skip1 = atoi( dosage_opt->getValue("skip1").c_str() );
  if ( dosage_opt->getValue("skip2") != "" ) 
    skip2 = atoi( dosage_opt->getValue("skip2").c_str() );

  if ( dosage_opt->isSet("dose1") )
    dosageScale2 = false;

  // Get relevant columns codes

  int snp_field = skip0 + 0;
  int a1_field = snp_field + skip1 + 1;
  int a2_field = a1_field + 1;
  int geno_field = a2_field + skip2 + 1;
  int pre_fields = 3 + skip0 + skip1 + skip2;
  
  if ( par::dosage_hasMap )
    {
      printLOG("A MAP file has been specified: only these markers will be processed\n");

      if ( par::dosage_hard_call )
	{
	  printLOG("Going to make hard-calls at " 
		   + dbl2str( par::dosage_hard_call_thresh ) + " threshold\n");
	  printLOG("and write remaining dosages to [ " + par::output_file_name 
		   + ".dosage" );
	  if ( compressed_out ) 
	    printLOG(".gz ]");
	  else
	    printLOG(" ]");
	  printLOG(" if more than " + int2str(par::dosage_hard_call_thresh2) + " non-calls\n");
	  par::SNP_major = true;


	  // Consider each new variant
	  // Space will already have been added
	  for (int l=0; l<nl_all; l++)
	    {	
	      // Causal variant, data for all individuals
	      // default to missing genotype
	      CSNP * newlocus = SNP[l];
	      newlocus->one.resize(n,true);
	      newlocus->two.resize(n,false);	      
	    }
	}
    }


  double thresh1_AA = (1-par::dosage_hard_call_thresh)/2.0;
  double thresh1_AB = 0.5 - ((1-par::dosage_hard_call_thresh)/2.0);
  double thresh2_AB = 0.5 + ((1-par::dosage_hard_call_thresh)/2.0);
  double thresh1_BB = 1 - (1-par::dosage_hard_call_thresh)/2.0;

    
  ////////////////////////////////////
  // Genotype dosage format options

  bool oneDose = false;
  bool twoProbs = false;
  bool threeProbs = false;

  if ( dosage_opt->getValue("format") == "1" ) 
    oneDose = true;
  else if ( dosage_opt->getValue("format") == "3" ) 
    threeProbs = true;
  else
    twoProbs = true;

  if ( oneDose )
    printLOG("Format set to one dosage per genotype\n");
  else if ( twoProbs )
    printLOG("Format set to two genotype probabilities\n");
  else
    printLOG("Format set to three genotype probabilities\n");

  int step = 1;
  if ( twoProbs ) step = 2;
  else if ( threeProbs ) step = 3;
   
  
  // Still assume a single FAM file; but allow for people/SNPs to be
  // spread across multiple files

  
  // *** Currently, the order still needs to be similar across files in the
  //     same SNP batch

  vector<string> dosageFilename_all;
  vector<string> headerFilename_all;
  vector<int> batchName;
  set<int>    batchNameSet;

  if ( ! filelist ) 
    {
      dosageFilename_all.push_back( par::dosage_file );
    }
  else
    {
      ZInput IN1( par::dosage_file , false );

      bool batchModeSet = false;
      
      while ( ! IN1.endOfFile() )
	{
	  vector<string> f = IN1.tokenizeLine();
	  if ( f.size() == 0 ) break;
	  
	  // possible formats, if !sepHeader
	  //    dosage-file
	  //    snp-batch dosage-file 
	  
	  // if sepHeader
	  //    dosage-file header-file
	  //    snp-batch dosage-file header-file
	  
	  
	  if ( batchModeSet ) 
	    {
	      if ( sepHeader && snpBatch && f.size() != 3 )
		error("Expecting 3 entries: batch dosage header");
	      else if ( sepHeader && (!snpBatch) && f.size() != 2 )
		error("Expecting 2 entries: dosage header");
	      else if ( (!sepHeader) && snpBatch && f.size() != 2 )
		error("Expecting 2 entries: batch dosage header");
	      else if ( (!sepHeader) && (!snpBatch) && f.size() != 1 )
		error("Expecting 1 entry: dosage");	      
	    }
	  else
	    {
	      if (sepHeader)
		{ 
		  if (f.size() == 3 ) 
		    snpBatch = true;
		  else if (f.size() == 2) 
		    snpBatch = false;
		  else
		    error("Problem with dosage file list format");
		}
	      else
		{
		  if ( f.size() == 2 )
		    snpBatch = true;
		  else if ( f.size() == 1 )
		    snpBatch = false;
		  else 
		    error("Problem with dosage file list format");
		}

	      batchModeSet = true;
	    }


	  // Store information

	  int term = 0;

	  if ( snpBatch )
	    {
	      int nm;

	      if ( ! from_string<int>( nm, f[term] , std::dec ) )
		error("Problem reading SNP batch number");
	      
	      batchName.push_back( nm );
	      batchNameSet.insert( nm );

	      ++term;		
	    }

	  dosageFilename_all.push_back(f[term++]);
	  
	  if ( sepHeader )
	    headerFilename_all.push_back(f[term++]);

	}

      IN1.close();      

      printLOG("Expecting " 
	       + int2str( dosageFilename_all.size() ) 
	       + " total files, in " + int2str(batchNameSet.size()) 
	       + " distinct batches of SNPs\n"); 
    }
  

  
  ///////////////////////////////////////////
  //
  // Set up some basic things, headers, etc
  //
  ///////////////////////////////////////////
  
  map<string,int> msample;
  for (int i=0; i<n; i++)
    msample.insert(make_pair(sample[i]->fid+"_"+sample[i]->iid,i));

  map<string,int> mlocus;
  if ( par::dosage_hasMap )
    for (int l=0; l<nl_all; l++)
      mlocus.insert(make_pair(locus[l]->name,l));

  
  string ext = countOccur 
    ? ".occur.dosage" : par::write_dosage 
    ? ".out.dosage" : ".assoc.dosage";

  map<string,int> occur;

  ZOutput detout;
  if ( par::dosage_hard_call ) 
    {
      ext = ".dosage";
      detout.open( par::output_file_name + ".dosage.det" , false );
      std::ostringstream s2( std::stringstream::out );
      s2 << sw("CHR",4)
	 << sw("SNP",par::pp_maxsnp)
	 << sw("BP", 12)
	 << sw("A1", 4)
	 << sw("A2", 4)
	 << sw("MAF", 8)
	 << sw("INFO", 8)
	 << sw("ABOVE", 8)
	 << sw("BELOW", 8)
	 << sw("RATE", 8) << "\n";

      detout.write( s2.str() );

    printLOG("Writing additional information to [ " 
	       + par::output_file_name + ".dosage.det ]\n");
    }
  
  string tail = compressed_out ? ".gz" : "";
  ZOutput zout( par::output_file_name + ext + tail , compressed_out );
  
  if ( ! par::dosage_hard_call )
    printLOG("Writing results to [ " + par::output_file_name + ext + tail + " ]\n");


  if ( par::dosage_hard_call || par::write_dosage )
    {
      // Write header to dosage output file
      std::ostringstream s2( std::stringstream::out );
      s2 << "SNP A1 A2 ";
      for (int i=0; i<n; i++)
	if ( ! sample[i]->missing ) 		  
	  s2 << sample[i]->fid << " "
	     << sample[i]->iid << " ";
      s2 << "\n";
      zout.write( s2.str() );
    }


  // Do we need a header?

  if ( ! ( countOccur || par::dosage_hard_call || par::write_dosage ) )
    {
      string es = par::bt ? "OR" : "BETA";
      if ( par::dosage_hasMap )
	zout << sw("CHR",4) 
	     << sw("SNP",12)
	     << sw("BP",12)
	     << sw("A1",4) 
	     << sw("A2",4) 
	     << sw("FRQ",8) 
	     << sw("INFO",8)
	     << sw(es,8)
	     << sw("SE",8)
	     << sw("P",8) << "\n";
      else
	zout << sw("SNP" ,12)
	     << sw("A1" ,4) 
	     << sw("A2" ,4) 
	     << sw("FRQ" ,8) 
	     << sw("INFO" ,8) 
	     << sw(es ,8)
	     << sw("SE",8)
	     << sw("P" ,8) << "\n";
      
    }



  ///////////////////////////////////////////////////
  // Set up association model

  bool OLD_assoc_glm_without_main_snp = par::assoc_glm_without_main_snp;
  bool OLD_clist = par::clist;

  par::assoc_glm_without_main_snp = true;
  par::clist = true;
  
  // Add an extra covariate slot

  ++par::clist_number;
  clistname.resize( par::clist_number );
  clistname[ par::clist_number - 1 ] = "DOSAGE";
  for (int i=0; i<n; i++)
    sample[i]->clist.resize( par::clist_number );

  int term = par::clist_number -1;
  int vcount = 0;
  
  
  ///////////////////////////////////////////////////
  // Set up scoring procedure

  if ( par::score_risk )
    {
      if ( par::score_risk_on_qrange )
	sizeMatrix(scores,n,qthresh.size());
      else
	sizeMatrix(scores,n,1);
    }




  ///////////////////////////////////////////
  //
  // Start looping through SNP batches
  //
  ///////////////////////////////////////////
  
  set<int>::iterator bi = batchNameSet.begin();
    
  while (1)
    {


      // Pull out the relevant set of files
      
      vector<string> dosageFilename;
      vector<string> headerFilename;
      
      if ( snpBatch )
	{
	  for (int i=0; i<dosageFilename_all.size(); i++)
	    {

	      if ( batchName[i] == *bi )
		{
		  dosageFilename.push_back( dosageFilename_all[i] );
		  
		  if ( sepHeader )
		    headerFilename.push_back( headerFilename_all[i] );
		  
		}
	    }
	}
      else
	{
	  dosageFilename = dosageFilename_all;
	  headerFilename =  headerFilename_all;
	}
      
            
      int nFiles = dosageFilename.size();
      vector<int> expected( nFiles );
      vector<ZInput*> vzin( nFiles );
      vector<ifstream*> vhead( nFiles );
      vector< vector<Individual*> > personMap( nFiles);
      vector<int> npeople( nFiles );

      int found = 0;
      int totpeople = 0;
      set<Individual*> inDosage;
      
      if ( nFiles > 1 && ! header ) 
	error("Can only specify noheader when reading a single dosage file\n");
      

      for (int f = 0 ; f < vzin.size() ; f++ )
	{
	  
	  ZInput * d = new ZInput( dosageFilename[f] , compressed_in );
	  vzin[f] = d;
	  
	  if ( sepHeader )
	    vhead[f] = new ifstream(headerFilename[f].c_str() , ios::in );
	  
	  // Read header, or use FAM file
	  
	  if ( header )
	    {
	      
	      vector<string> tok;
	      
	      if ( ! sepHeader )
		tok = vzin[f]->tokenizeLine();
	      else
		{
		  // read all ID pairs
		  while (1)
		    {
		      string fid, iid;
		      (*vhead[f]) >> fid >> iid;
		      if ( fid == "" || vhead[f]->eof() )
			break;
		      tok.push_back(fid);
		      tok.push_back(iid);		  
		    }	     
		}
	      
	  int firstCol = sepHeader ? 0 : pre_fields ;
	  if ( tok.size() < firstCol ) 
	    error("Bad format fdr dosage file, expecting more columns");
	  if ( ! sepHeader )
	    {
	      if ( tok[snp_field] != "SNP" || tok[a1_field] != "A1" || tok[a2_field] != "A2" )
		error("Badly aligned columns for: SNP A1 A2");
	    }
	  
	  if ( (tok.size() - firstCol ) % 2 != 0 )
	    error("Expecting 3 + 2 * N columns in header\n");

	  // Based on header field, two entries per person (FID, IID)
	  npeople[f] = (tok.size()- firstCol )/2;
	  	 
	  expected[f] = pre_fields + npeople[f] * step;

	  
	  //////////////////////////////////////////////////////
	  // Process header row

	  for (int i=firstCol; i<tok.size(); i+= 2 )
	    {
	      	      
	      string id = tok[i] + "_" + tok[i+1];
	      
	      map<string,int>::iterator m = msample.find(id);
	      
	      if ( m == msample.end() )
		{
		  personMap[f].push_back( (Individual*)NULL ) ;
		}
	      else
		{
		  if ( inDosage.find( sample[ m->second ] ) != inDosage.end() )
		    error("The person appears in >1 dosage file: " + id );
		  Individual * person = sample[m->second];
		  personMap[f].push_back( person );
		  inDosage.insert( person );
		  ++found;
		}
	    }
	  
	  totpeople += npeople[f];
	}
      else
	{
	  // If no explicit header given
	  npeople[f] = n;
	  expected[f] = pre_fields + npeople[f] * step;
	  totpeople += npeople[f];
	  for (int i=0; i<n; i++ )
	    {
	      string id = sample[i]->fid+ "_" + sample[i]->iid;
	      map<string,int>::iterator m = msample.find(id);
	      if ( m == msample.end() )
		{
		  personMap[f].push_back( (Individual*)NULL ) ;
		}
	      else
		{
		  if ( inDosage.find( sample[ m->second ] ) != inDosage.end() )
		    error("The person appears in >1 dosage file: " + id );
		  Individual * person = sample[m->second];
		  personMap[f].push_back( person );
		  inDosage.insert( person );
		}
	    }
	}
	  
	  
	} // next dosage file
      

      if ( ! snpBatch  && header ) 
	{
	  printLOG("Matched to " + int2str(found) 
		   + " of " + int2str( totpeople ) + " individuals");
	  if ( filelist ) 
	    printLOG(" in " + int2str( nFiles ) + " files\n");
	  else
	    printLOG(" in [ " + par::dosage_file + " ]\n");
	}
      
      
      if ( sepHeader )
	{
	  for (int f = 0 ; f < nFiles; f++ ) 
	    vhead[f]->close();
	}
      
      // Remove missing individuals
      
      if ( ! par::dosage_hard_call )
	{
	  int n_removed = keepIndividuals( inDosage );
	  
	  if ( n_removed > 0 )
	    {
	      printLOG("Removed " + int2str(n_removed) 
		       + " individuals not in dosage file\n");	  
	      // Update person map, given we've removed some people
	      msample.clear();
	      for (int i=0; i<n; i++)
		msample.insert(make_pair(sample[i]->fid+"_"+sample[i]->iid,i));	  
	    }   
	}
      
      
      // Create final file column position -> sample # mapping, for each file
      
      vector< vector<int> > personPosition( nFiles );
      for (int f = 0 ; f < nFiles; f++ )
	{
	  
	  for (int i = 0; i < personMap[f].size(); i++)
	    {
	      Individual * person = personMap[f][i]; 
	      if ( person == NULL )
		{
		  personPosition[f].push_back(-1);
		}
	      else
		{
		  string id = person->fid + "_" + person->iid;
		  map<string,int>::iterator p = msample.find( id );
		  personPosition[f].push_back( p->second );
		}
	    }
	}
      


      /////////////////////////////////////////////////
      // Read main dosage data, and analyse SNP by SNP
      
      while ( 1 ) 
	{
	  
	  bool done = false;
	  bool skip = false;
	  
	  string snp_id;
	  string a1_id; 
	  string a2_id;
	  int snp_code;
	  
	  int goodCall = 0;
	  int badCall = 0;
	  map<Individual*,double> dose1;
	  map<Individual*,double> dose2;
	  
	  for ( int f = 0 ; f < nFiles ; f++ )
	    {
	      
	      if ( vzin[f]->endOfFile() )
		{
		  done = true;
		  break;
		}
	      
	      // Read line from the correct file
	      vector<string> tok = vzin[f]->tokenizeLine();
	      
	      if ( tok.size() == 0 )
		{
		  skip = true; 
		  break;
		}
	      
	      if ( tok.size() != expected[f] )
		error("Problem with line:\n" + displayLine(tok) );
	      
	      
	      if ( tok[ snp_field ] == "" ) 
		{
		  skip = true;
		  break;
		}
	      
	      if ( f == 0 ) 
		{
		  snp_id = tok[ snp_field ];
		  a1_id = tok[ a1_field ];
		  a2_id = tok[ a2_field ];
		}
	      else
		{
		  if ( snp_id != tok[ snp_field ] )
		    error("Misaligned SNPs in dosage files");
		  if ( a1_id != tok[a1_field ] ) 
		    error("Misaligned allele codes in dosage file");
		}
	      
	      
	      // If we've loaded a MAP file, then ignore this
	      // marker if it is not present in the MAP file
	      
	      if ( par::dosage_hasMap )
		{
		  map<string,int>::iterator mi = mlocus.find( snp_id );
		  if ( mlocus.find( snp_id ) == mlocus.end() )
		    {
		      skip = true;
		      continue;		  
		    }
		  
		  if ( f==0 )
		    {
		      snp_code = mi->second;
		      locus[ snp_code ]->allele1 = a1_id;
		      locus[ snp_code ]->allele2 = a2_id;
		    }
		}
	      
	      
	      // Are we just in the mode in which we simply count
	      // how many times/files we see this SNP?
	  
	      if ( countOccur ) 
		{
		  map<string,int>::iterator o = occur.find( tok[ snp_field ] );
		  if ( o != occur.end() ) 
		    (o->second)++;
		  else
		    occur.insert(make_pair( tok[ snp_field ], 1 ));
		  continue;
		}
	      
	      // Start position for genotype data
	      
	      int j = pre_fields;
	      
	      for (int i=0; i<npeople[f]; i++)
		{
		  
		  bool problem = false;
		  
		  // Skip this individual if not in file
		  if ( personPosition[f][i] == -1 )
		    {
		      j += step;
		      continue;
		    }


		  Individual * person = sample[ personPosition[f][i] ];

		  double dose, d1, d2, d3;
		  
		  if ( oneDose ) 
		    {
		      if ( !from_string<double>( dose, tok[j++], std::dec ) )
			problem=true;
		      if ( dosageScale2 )
			{
			  if ( dose < 0 || dose > 2 ) 
			    problem = true;	      
			  dose /= 2.0;
			}
		      else
			{
			  if ( dose < 0 || dose > 1 ) 
			    problem = true;	      		      
			}
		    }
		  else if ( twoProbs )
		    {
		      if ( !from_string<double>( d1, tok[j++], std::dec ) )
			problem = true;
		      if ( !from_string<double>( d2, tok[j++], std::dec ) )
			problem = true;
		      if ( d1 < 0 || d1 > 1 ) 
			problem = true;
		      if ( d2 < 0 || d2 > 1 ) 
			problem = true;
		      if ( d1 + d2 > 1 ) 
			problem = true;	   
		      dose = d1 + d2/2.0;
		      dose1.insert(make_pair( person , d1 ) );
		      dose2.insert(make_pair( person , d2 ) );
		    }
		  else if ( threeProbs )
		    {
		      if ( !from_string<double>( d1, tok[j++], std::dec ) )
			problem = true;
		      if ( !from_string<double>( d2, tok[j++], std::dec ) )
			problem = true;
		      if ( !from_string<double>( d3, tok[j++], std::dec ) )
			problem = true;
		      if ( d1 < 0 || d1 > 1 ) 
			problem = true;
		      if ( d2 < 0 || d2 > 1 ) 
			problem = true;	      
		      // skip sanity check on 3rd
		      // 	      if ( d1 + d2 > 1 ) 
		      // 		problem = true;	   
		      dose = d1 + d2/2.0;
		      
		      dose1.insert(make_pair( person , d1 ) );
		      dose2.insert(make_pair( person , d2 ) );
		      
		    }
		  
		  // Do we want to make a hard call now, or store 
		  // dosage for analysis?
		  if ( par::dosage_hard_call )
		    {
		      bool s1 = true;
		      bool s2 = false;
		      
		      if ( oneDose )
			{
			  if ( dose < thresh1_AA ) 
			    { s1=s2=false; }
			  else if ( dose > thresh1_AB && dose < thresh2_AB )
			    { s1=false; s2=true; }
			  else if ( dose > thresh1_BB ) 
			    { s1=s2=true; }
			}
		      else
			{
			  if ( d1 > par::dosage_hard_call_thresh ) 
			    { s1=s2=false; }
			  else if ( d2 > par::dosage_hard_call_thresh )
			    { s1=false; s2=true; }
			  else if ( 1-d1-d2 > par::dosage_hard_call_thresh )
			    { s1=s2=true; }
			}
		      
		      if ( s1 && ! s2 ) 
			++badCall;
		      else
			++goodCall;
		      
		      SNP[snp_code]->one[personPosition[f][i]] = s1;
		      SNP[snp_code]->two[personPosition[f][i]] = s2;		      
		      
		    }
		  

		  if ( problem )
		    {
		      sample[ personPosition[f][i] ]->missing2 = true;
		    }
		  else
		    {
		      sample[ personPosition[f][i] ]->missing2 = false;
		      sample[ personPosition[f][i] ]->clist[ term ] = dose;
		      
		    }
		  
		}
	      
	    } // Next dosage file
	  

	  // Are we at the end of a file though?
	  if ( done ) 
	    break;
	  if ( skip ) 
	    continue;
	  
	  // Do we need to bother processing the genotype dosage data?
	  if ( countOccur ) 
	    continue;
	  
	  ++vcount;
	  if ( ! par::silent ) 
	    cerr << "Processed " << vcount << " markers                     \r";
	  
      

	  
	  
	  ///////////////////////////////////////////////////
	  // Set up scoring procedure
      
	  if ( par::score_risk )
	    {
	  
	      // Does this variant have a score?
	      
	      Var v;
	      v.snp = snp_id;
	      map<Var,double>::iterator i = wt.find( v );
	      if ( i == wt.end() )
		continue;
	      double weight = i->second;
	      
	      // Right allele?
	      bool swapAllele = false;	 
	      if ( i->first.a1 == a2_id )
		{
		  // need to swap
		  swapAllele = true;	      
		}
	      else
		{
		  if ( i->first.a1 != a1_id )
		    continue;
		}
	      
	      if ( par::score_risk_on_qrange )
		{
		  map<Var,vector<int> >::iterator i0 = qflag.find(v);
		  if ( i0 == qflag.end() )
		    continue;
		  vector<int> & inQ = i0->second;
		  for (int q=0; q<inQ.size(); q++)
		    {
		      for (int i=0; i<n; i++)
			scores[i][ inQ[q] ] += swapAllele ? 
			  weight * ( 1 - sample[i]->clist[ term ] ) :
			  weight * sample[i]->clist[ term ];
		    }
		}
	      else
		{
		  for (int i=0; i<n; i++)
		    {
		      scores[i][0] += swapAllele ? 
			weight * ( 1 - sample[i]->clist[ term ] ) : 
			weight * sample[i]->clist[ term ];
		    }
		}
	      
	      // Skip association test, etc
	      continue;
	    }
	  
	  
	  ///////////////////////////////////////////
	  // tabulate frequency, and info/r^2 score
	  
	  double frq = 0;
	  int cnt = 0;
	  for (int i=0; i<n; i++)
	    if ( ! sample[i]->missing ) 
	      {
		frq += sample[i]->clist[ term ];
		++cnt;
	      }
	  frq /= (double)cnt;
	  
	  double theoreticalVariance = frq * ( 1 - frq );
	  double dosageSSQ = 0;
	  for (int i=0; i<n; i++)
	    if ( ! sample[i]->missing ) 
	      {
		double t1 = sample[i]->clist[ term ] - frq ;
		t1 *= t1;
		dosageSSQ += t1;
	      }
	  
	  
	  double empiricalVariance  = 2 * ( dosageSSQ / (double)cnt);
	  double rsq = theoreticalVariance > 0 ? empiricalVariance / theoreticalVariance : 0;

	  

	  //////////////////////////////////////////
	  // Give some output

	  if ( par::dosage_hard_call )
	    {
	      std::ostringstream s2( std::stringstream::out );
	      s2 << sw(locus[snp_code]->chr, 4)
		 << sw(snp_id, par::pp_maxsnp) 
		 << sw(locus[snp_code]->bp, 12) 
		 << sw(a1_id, 4)   
		 << sw(a2_id, 4)   
		 << sw(frq,4, 8)   
		 << sw(rsq,4,8)  
		 << sw(goodCall, 8)  
		 << sw(badCall, 8)   
		 << sw((double)goodCall/(double)(goodCall+badCall),4, 8) << "\n";
	      
	      detout.write( s2.str() );

	      // Do we need to write this line back out as dosage info? 
	      if ( badCall > par::dosage_hard_call_thresh2 )   
		{

		  // Write header to dosage output file

		  std::ostringstream s2( std::stringstream::out );
		  
		  s2 << snp_id << " " << a1_id << " " << a2_id << " ";
		  
		  if ( oneDose )
		    {
		      for (int i=0; i<n; i++)	
			if ( ! sample[i]->missing ) 		  
			  s2 << 2 * sample[i]->clist[term] << " ";
		    }
		  else
		    {
		      for (int i=0; i<n; i++)	
			if ( ! sample[i]->missing )
			  {
			    double d1 = dose1.find(sample[i])->second;
			    double d2 = dose2.find(sample[i])->second;
			    if ( twoProbs )
			      s2 << d1 << " " << d2 << " ";
			    if ( threeProbs )
			      s2 << 1 - d1 - d2 << " ";
			  }
		    }
		  s2 << "\n";
		  zout.write( s2.str() );
		}
	      

	      continue;
	    }

	  
	  if ( par::write_dosage )
	    {

	      std::ostringstream s2( std::stringstream::out );
	      
	      s2 << snp_id << " " << a1_id << " " << a2_id << " ";
	      
	      if ( oneDose )
		{
		  for (int i=0; i<n; i++)	
		    if ( ! sample[i]->missing ) 		  
		      s2 << 2 * sample[i]->clist[term] << " ";
		}
	      else
		{
		  for (int i=0; i<n; i++)	
		    if ( ! sample[i]->missing )
		      {
			double d1 = dose1.find(sample[i])->second;
			double d2 = dose2.find(sample[i])->second;
			if ( twoProbs )
			  s2 << d1 << " " << d2 << " ";
			if ( threeProbs )
			  s2 << 1 - d1 - d2 << " ";
		      }
		}
	      
	      s2 << "\n";
	      
	      zout.write( s2.str() );
	      

	      // Do not perform association test, go straight to 
	      // next marker
	      
	      continue;
	      
	    }
	  

	  

	  ///////////////////////////////////////////
	  // Perform association
	  
	  glmAssoc(false,*pperm);
	  
      
	  ///////////////////////////////////////////
	  // Report results

	  bool valid = model->isValid();
      
	  // Do not output for bad markers
	  // Hard-code for now...
	  
	  if ( frq < 0.01 || frq > 0.99 || rsq < 0.1 || rsq > 2) 
	    {
	      valid = false;
	      if ( rsq<0 ) rsq = 0;
	      if ( rsq>2 ) rsq = 2.0;	  
	    }
	  
	  vector_t b = model->getCoefs();
	  vector_t pval = model->getPVals();
	  vector_t var = model->getVar();
	  
	  // NOTE: b includes intercept; pval doesn't
	  
	  // Note: internal coding of dosage is on 0..1 scale, so
	  // divide beta and SE by 2 here to get per-allele effects

	  double statistic = valid ? model->getStatistic() : 0;
	  double pvalue = pval[ pval.size()-1 ];
	  double beta = par::bt ? exp( b[ b.size()-1 ] / 2.0) : b[ b.size()-1 ] / 2.0 ;
	  
	  double se = sqrt( var[ var.size()-1 ] ) / 2.0;


	  
	  if ( par::dosage_hasMap )
	    zout << sw(locus[snp_code]->chr , 4)
		 << sw(snp_id , 12)
		 << sw(locus[snp_code]->bp ,12)
		 << sw(a1_id ,4)
		 << sw(a2_id ,4)
		 << sw(frq,4 ,8)
		 << sw(rsq,4 ,8);
	  else
	    zout << sw(snp_id ,12)
		 << sw(a1_id ,4)
		 << sw(a2_id ,4)
		 << sw(frq,4 ,8)
		 << sw(rsq,4 ,8);
	  
	  if ( valid ) 
	    {
	      zout << sw(beta,4,8)
		   << sw(se,4,8)
		   << sw(pvalue,-4,8)<< "\n";
	    }
	  else
	    {
	      zout << sw("NA",8) 
		   << sw("NA",8)
		   << sw("NA",8) << "\n";
	    }
	  
      
	  delete model;
	  
	  // Next variant(s)
	}
      
        
      ///////////////////////////////////////
      //
      // Finished this batch
      //
      ///////////////////////////////////////
      
      for (int f = 0 ; f < nFiles; f++)
	vzin[f]->close();
      
      if ( ! snpBatch ) 
	break;
      
      ++bi;
      
      if ( bi == batchNameSet.end() )
	break;
      
    }
  

  ///////////////////////////////////////
  //
  // Finished all .. output & wrap up
  //
  ///////////////////////////////////////

  if ( ! par::silent ) 
    cerr << "\n";

  // Write out occurence info
  
  if ( countOccur ) 
    {
      map<string,int>::iterator o = occur.begin();
      int totCount = 0;
      int nonF = 0;
      while ( o != occur.end() )
	{
	  zout.write( o->first + " " + int2str( o->second ) + "\n" );
	  if ( o->second != dosageFilename_all.size() )
	    ++nonF;
	  totCount += o->second;
	  ++o;
	}
      printLOG("Counted unique " + int2str(occur.size() ) 
	       + " markers, " + int2str( totCount ) 
	       + " across all files\n");
      if ( nonF == 0 ) 
	printLOG("All SNPs occured exactly once in all files\n");
      else
	printLOG("Note: " + int2str(nonF) + " SNPs did not occur exactly once per file\n");
      
    }

  
  zout.close();
  
  if ( par::dosage_hard_call ) 
    detout.close();

//    if ( ! countOccur )
//      printLOG("In total, analysed " + int2str(vcount) + " markers\n");



  ///////////////////////////////////////
  // Write scores to file

  if ( par::score_risk )
    {
      int qq = 0;
      while (1)
	{
	  string append = par::score_risk_on_qrange ? 
	    ".S" + int2str(qq+1) : "";

	  ofstream O1( ( par::output_file_name + append + ".profile").c_str() , ios::out );	  

	  O1 << setw(par::pp_maxfid) << "FID" << " " 
	     << setw(par::pp_maxiid) << "IID" << " "
	     << setw(6) << "PHENO" << " " 
	     << setw(8) << "SCORE" << "\n";
	  
	  for ( int i=0; i<n; i++ )
	    {
	      Individual * person = sample[i];	      
	      O1 << setw(par::pp_maxfid) << person->fid << " " 
		 << setw(par::pp_maxiid) << person->iid << " "
		 << setw(6) << person->phenotype << " " 
		 << setw(8) << scores[i][qq] << "\n";
	    }
  
	  O1.close();

	  if ( !par::score_risk_on_qrange )
	    break;

	  if ( ++qq == qthresh.size() )
	    break;
	}
    }



  ///////////////////////////////////////
  // Some final tidying up
  
  par::assoc_glm_without_main_snp = OLD_assoc_glm_without_main_snp;
  par::clist = OLD_clist;  
  --par::clist_number;
  clistname.resize( par::clist_number );
  for (int i=0; i<n; i++)
    sample[i]->clist.resize( par::clist_number );
  
      
  return;
  
}




void setUpQScoring( map<string,double> & qscore,
		    vector<double2> & qthresh,
		    vector<string> & qlabel)
{
      
  checkFileExists( par::score_qfile );
  checkFileExists( par::score_qrange_file );
  
  PP->printLOG("Reading quantitative scores from [ " + par::score_qfile + " ]\n");
  PP->printLOG("Reading score ranges from [ " + par::score_qrange_file + " ]\n");
  
  ifstream Q1( par::score_qfile.c_str() , ios::in );
  while ( ! Q1.eof() )
    {
      string snp;
      string str_score;
      double score;	  
      Q1 >> snp >> str_score;	  
      if ( ! from_string<double>( score , str_score , std::dec ) )
	continue;
      if ( snp == "" )
	continue;
      
      qscore.insert( make_pair( snp , score ) );
    }

  Q1.close();
  PP->printLOG("Read q-scores for " + int2str( qscore.size() ) + " SNPs\n");      
  
  Q1.open( par::score_qrange_file.c_str() , ios::in );
  while ( ! Q1.eof() )
    {
      // Expect: name, lower, upper
      string label;
      double lower, upper;
      Q1 >> label >> lower >> upper;
      if ( label == "" )
	continue;
      double2 d2(lower,upper);
      qthresh.push_back( d2 );
      qlabel.push_back( label );
    }
  Q1.close();
  PP->printLOG("Read " + int2str( qthresh.size() ) + " thresholds to apply\n");            
}



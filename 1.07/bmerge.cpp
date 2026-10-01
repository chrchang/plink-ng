

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2008 Shaun Purcell                  //
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
#include <map>
#include <algorithm>
#include <bitset>
#include "plink.h"
#include "options.h"
#include "helper.h"


void Plink::mergeBinaryData()
{
  
  // Function to merge a text file with an exisiting data set
  // either SNP-major or individual-major modes
  
  if (!par::merge_list)
    {
      printLOG( "Using merge mode " + int2str( par::merge_mode ) + " : ");
      if (par::merge_mode==1) printLOG("consensus call (default)\n");
      else if (par::merge_mode==2) printLOG("overwrite if missing in original\n");
      else if (par::merge_mode==3) printLOG("overwrite unless missing in new\n");
      else if (par::merge_mode==4) printLOG("overwrite none\n"); 
      else if (par::merge_mode==5) printLOG("overwrite all\n");
      else if (par::merge_mode==6) printLOG("diff mode: all differences\n");
      else if (par::merge_mode==7) printLOG("diff mode: non-missing differences\n");

      diff_overlap = 0;
      diff_nonmissing_overlap = 0;
      diff_concordant_overlap = 0;

    }
  
  // We've already loaded in the first file
  // Do not overwrite any existing phenotype information
  
  checkFileExists(par::merge_bedfile);
  checkFileExists(par::merge_bimfile);
  checkFileExists(par::merge_famfile);
  
  // Make hash of original SNP names
  
  map<string,int> mlocus;
  for (int l=0;l<nl_all;l++)
    mlocus.insert(make_pair(locus[l]->name,l));
  map<string,int>::iterator ilocus;
  
  
  // A temporary hash for the names of any markers that
  // do not match in terms of strand
  set<string> misstrand;
  map<string,int> misstrand_dummy;
  bool fatal_error = false;


  ///////////////////////////////////////
  // .bim 

  vector<Locus> ordered(0);
  map<string,bool> exists;
  vector<Locus*> locus2(0);
  set<string> flip_alleles;

  ifstream MAP(par::merge_bimfile.c_str(), ios::in);
  MAP.clear();

  int exist_cnt=0;
  
  int c=0;
  while(!MAP.eof())
    {
      
      Locus * loc = new Locus;

      long int inc;

      MAP >> loc->chr   // will automatically by numeric
	  >> loc->name 
	  >> loc->pos   // will automatically be in M units
	  >> loc->bp     
	  >> loc->allele1
	  >> loc->allele2;
      
                  
      inc = loc->bp;
            
      // Use the frequency slot temporarily to 
      // store order information
      loc->freq = c++;
      
      // Check that cM/M specification looks correct, if 
      // we want to perform a plink-based analysis
      if (par::plink && (!par::cm_map) && (loc->pos > 50) )
	error("Looks like you need to specify --cm ??");
      
      // Convert cM to M map distances
      if (par::cm_map) loc->pos /= 100;

      // Including all loci in merge-mode
      if (loc->name!="") 
	{
	  
	  ilocus = mlocus.find(loc->name);
	      

	  ///////////////////////////////////////////////////
	  // Check whether or not this Locus already exists?
	  
	  if (ilocus != mlocus.end())
	    {
	      
	      Locus * loc2 = locus[ilocus->second];
	      
	      // Check same chromosome and positions, etc
	      if ( loc2->chr != loc->chr )
		{
		  cerr << "Warning: different chromosome for " 
		       << loc->name << "\n";
		  loc->chr = loc2->chr;
		}
	      
	      // Check same chromosome and positions, etc
	      if ( loc2->bp != loc->bp )
		{
		  cerr << "Warning: different physical position for " 
		       << loc->name << "\n";
		  loc->bp = loc2->bp;		  
		}
	      
	      if ( loc2->pos != loc->pos )
		{
		  cerr << "Warning: different genetic position for " 
		       << loc->name << "\n";
		  loc->pos = loc2->pos;		  
		}
	      
	      exists[loc->name] = true;
	      exist_cnt++;
	      locus2.push_back(loc2);
	      
	      // Keep the new file order (would have been in freq)
	      int t = (int)loc->freq;
	      
	      if ( loc2->allele1 == "" ) 
		loc2->allele1 = par::missing_genotype;
	      
	      if ( loc2->allele2 == "" ) 
		loc2->allele2 = par::missing_genotype;
	      
	      /////////////////////////////////////////
	      // Add allele names to list, if needed
	      // and check if allele codes need flipping?
	      // e.g. if A/C SNP in memory but C/A in file?
	      // or had one allele missing
	      // and check that alleles match up correctly
	      
	      // 0 A   A 0  ->   0 A + flip 
	      // 0 0   0 0  ->   0 0 
	      // 0 0   0 A  ->  {0 A}
	      // 0 0   A B  ->  {A B}
	      // 0 A   0 A  ->   0 A
	      // 0 A   0 B  ->  {B A} + flip
	      // 0 A   A B  ->  {B A} + flip
	      // 0 A   B A  ->  {B A} 
	      // 0 A   0 0  ->   0 A
	      // A B   A B  ->   A B
	      // A B   B A  ->   A B + flip
	      // A B   0 A  ->   A B + flip 
	      // A B   0 0  ->   A B
	      
	      // 1) Are there any empty slots in the existing allele? 
	      //    If so, fill them in with any new alleles
	      
	      // 2) Do we have a strand problem? 
	      
	      // 3) Or do we need a flip? 
	      
	      
	      // New codes
	      
	      string one = loc->allele1;
	      string two = loc->allele2;
	            
//     	      cout << "LOCUS " << loc->name << " OLD, NEW = [" 
//     		   << loc2->allele1 << "] [" 
//     		   << loc2->allele2 << "]    [" 
//     		   << one << "] [" 
//     		   << two << "]\n"; 
	      
	      set<string> alleleCount;

	      if ( one != par::missing_genotype )
		alleleCount.insert( one );
	      if ( two != par::missing_genotype )
		alleleCount.insert( two );
	      if ( loc2->allele1 != par::missing_genotype )
		alleleCount.insert( loc2->allele1 );
	      if ( loc2->allele2  != par::missing_genotype )
		alleleCount.insert( loc2->allele2 );
	      
	      // More than 2 obseved alleles? 

	      if ( alleleCount.size() > 2 ) 
		{
		  misstrand.insert(loc2->name);
		  fatal_error = true;
		}
	      else
		{
		  
		  //////////////////////////////
		  // 1) Fill in empty slots
		  
		  // Fill slot 2 first (i.e. 0 A code for monomorphic, not A 0)
		  
		  // If first new allele is not missing...
		  if (one!=par::missing_genotype)
		    {
		      // ...and not already listed
		      if (one!=loc2->allele1 && one!=loc2->allele2)
			{
			  // ...then add to first empty slot
			  if(loc2->allele2=="" || loc2->allele2==par::missing_genotype)
			    loc2->allele2=one;
			  else if(loc2->allele1=="" || loc2->allele1==par::missing_genotype)
			    loc2->allele1=one;
			}
		    }
		  
		  if (two!=one && two!=par::missing_genotype)
		    // ...and not already listed
		    if (two!=loc2->allele1 && two!=loc2->allele2)
		      {
			// ...then add to first empty slot
			if(loc2->allele2=="" || loc2->allele2==par::missing_genotype) 
			  loc2->allele2=two;
			else if(loc2->allele1=="" || loc2->allele1==par::missing_genotype ) 
			  loc2->allele1=two;
		      }
		  
		  
		  //////////////////////////////
		  // 2) Need a flip?
		  
		  if ( ( one == loc2->allele2 && one != par::missing_genotype ) ||
		       ( two == loc2->allele1 && two != par::missing_genotype ) ||
		       ( one != loc2->allele1 && one != par::missing_genotype && loc2->allele1 != par::missing_genotype ) || 
		       ( two != loc2->allele2 && two != par::missing_genotype && loc2->allele2 != par::missing_genotype ) )
		    {
		      
		      if ( one == loc2->allele2 || two == loc2->allele1 ) 
			{
			  flip_alleles.insert(loc2->name);
			  
			}
		      else 
			{
			  
			  //////////////////////////////
			  // 3) Strand, wrong coding?
			  
			  misstrand.insert(loc2->name);
			  fatal_error = true;
			}
		    }
	      
		  
		}
	      
	      // Clean up what we do not need
	      delete loc;

	      // Replace with old locus (but swap back in new file position)
	      loc = loc2;
	      loc->freq = t;
	      
	    }
	  else 
	    {
	      // Locus does not exist -- add to locus list
	      exists[loc->name] = false;
	      locus2.push_back(loc);
	    }
	  
	  
	  ordered.push_back(*loc);
	}
      
      
    }

  MAP.clear();
  MAP.close();
  

   
  ///////////////////////////////////////////////////////
  // Did we encounter any fatal errors from flipped SNPs? 
   
  if (fatal_error)
    {
      
      ofstream MSNP;
      string f = par::output_file_name+".missnp";
      MSNP.open(f.c_str(), ios::out);
      set<string>::iterator ilocus;
      
      for ( ilocus = misstrand.begin() ; ilocus != misstrand.end() ; ilocus++)
	{
	  MSNP << *ilocus << "\n";
	}
      MSNP.close();
      printLOG("\nFound " + int2str(misstrand.size()) + " SNPs that do not match in terms of allele codes\n");
      printLOG("Might include strand flips, although flipped A/T and C/G SNPs will be undetected)\n");
      printLOG("Writing problem SNPs to [ " + f + " ]\n");
      error("Stopping due to mis-matching SNPs -- check +/- strand?");
    }
  

  
  if (!par::merge_list)
    {
      printLOG("\n" +int2str(locus2.size()) + 
	       " markers to be merged from [ " +par::merge_bimfile + " ]\n");
      printLOG("Of these, "+int2str( locus2.size()-exist_cnt ) + " are new, " +
	       int2str( exist_cnt ) +  " already exist in current data\n");
    }



  ///////////////////////////////////////////////
  // Build ordered table, so that genotypes can be inserted 
  // in correct order; then swap locus file over
  
  // Sorting a vector of pointers, so we need this special fix
  stable_sort(locus2.begin(),locus2.end(),less<Locus*>());
  
  // Sorting a normal vector
  stable_sort(ordered.begin(),ordered.end());
  
  c=0;
  for (int i=0; i<ordered.size(); i++)
    {
      // swap file order into locus position
      // and make all same chromosome
      ordered[i].bp = (int)ordered[i].freq;
      ordered[i].chr = 1;
      
	// keep track of genetic order, but only 
	// for nonmissing loci
	ordered[i].freq = c++;
    }
  
  // resort to get lookup table
  stable_sort(ordered.begin(),ordered.end());  

  // i.e. for file position k, the locus position is ordered[k]->freq

  // p2 p3 p1 p5 p4   : genetic position
  // 0  1  2  3  4    : file order

  // sort by cM
  // p1 p2 p3 p4 p5   : genetic 
  // 2  0  1  4  3    : file order
  // 0  1     2       : add genetic order: nonmissing...
  // 

  // sort by file order again
  // p2 p3 p1 p5 p4   : genetic
  // 0  1  2  3  4    : file
  // 1     0     2    : position to put in locus[l]
  

  // Add new locus2() to end of locus()
  for (int l=0; l<locus2.size(); l++)
    {
      // need to add a new MAP entry
      if (!exists.find(locus2[l]->name)->second )
	{
	  Locus * loc = new Locus;
	  loc = locus2[l];
	  locus.push_back(loc);
	}
    }
  

  
  ///////////////////////////////////////////////
  // .fam
  
  // Make new hash of Locus names
  
  mlocus.clear();
  for (int l=0;l<locus.size();l++)
    {
      mlocus.insert(make_pair(locus[l]->name,l));
    }
  if (mlocus.size() != locus.size() ) 
    {
      cerr << "Problem encountered merging files, with the following markers:\n";
       mlocus.clear();
       for (int l=0;l<locus.size();l++)
 	{
 	  if (mlocus.find(locus[l]->name) != mlocus.end())
 	    cerr << locus[l]->name << "\n";
 	  mlocus.insert(make_pair(locus[l]->name,l));
 	}
       cerr <<  "[ dump info: sizes = " << mlocus.size() << " and " << locus.size() << " ]\n";
       error("Cannot merge files. Check your MAP files.");
    } 
  
  // Make hash of existing individuals
  map<string,int> msample;
  for (int i=0;i<n; i++)
    msample.insert(make_pair(sample[i]->fid+"_"+sample[i]->iid,i));
  map<string,int>::iterator isample;


  // Resize all existing individuals
  // and set new elements to missing genotype (TF)
      
  if (par::SNP_major)
    {
      // Add space for new SNPs
      for (int i=0; i<locus2.size()-exist_cnt; i++)
	{
	  CSNP * newlocus = new CSNP;
	  newlocus->one.resize(n,true);
	  newlocus->two.resize(n,false);
	  SNP.push_back(newlocus);
	}      
      
    }
  else
    {
      // If using individual-major mode 
      for (int i=0; i<n; i++)
	{
	  sample[i]->one.resize(locus.size(),true);
	  sample[i]->two.resize(locus.size(),false);
	}
    }
  
  


  // An output file for diff mode
  ofstream MERD;
  if (par::merge_mode >=6)
    {
      string f = par::output_file_name+".diff";
      MERD.open(f.c_str(), ios::out);
      MERD << setw(20) << "SNP" << " " 
	   << setw(20) << "FID" << " " 
	   << setw(20) << "IID" << " " 
	   << setw(8) << "NEW" << " " 
	   << setw(8) << "OLD" << " " << "\n";
    }
  
  int new_person = 0;
  int old_person = 0;
  
  
  ///////////////////////////////////////
  // Read in FAM/BED file to new merge file
  

  // Initially, assume a binary trait
  par::qt = false;
  par::bt = true; 
  
  
  ifstream FAM;
  FAM.open(par::merge_famfile.c_str());
  FAM.clear();
  
  vector<bool> existing_person_list;
  int original_sample_size = sample.size();

  vector<Individual*> sample2;

  c=0;
  while(!FAM.eof())
    {

      Individual * person = new Individual;

      // No comments allowed in BED/BIM/FAM files

      // First 6 obligatory fields
      string phenotype;

      FAM >> person->fid 
	  >> person->iid 
	  >> person->pat
	  >> person->mat
	  >> person->sexcode
	  >> phenotype;
      

      // Skip last empty line that gets read
      if (person->fid=="") break;
      
      // Check for reserved family ID code
      if ( person->fid=="FID" )
	error("FID is a reserved ID... please select a different family ID");
      

      // Are we using 0/1 coding?
      if (par::coding01) 
	{
	  if ( phenotype == "1" ) 
	    phenotype = "2";      
	  else if ( phenotype == "0" )
	    phenotype = "1";
	  else 
	    phenotype = "0";
	}
            
      if (person->sexcode=="1")
	person->sex = true;  // male
      else if (person->sexcode=="2")
	person->sex = false; // female (default)
      else if (!par::ignore_missing_sex)
	{
	  person->missing = true;
	}

      
      // Have we already created this person?

      bool already_in = false;
      isample = msample.find(person->fid+"_"+person->iid);
      int indn = isample->second;
      if ( isample != msample.end() )
	{
	  already_in = true;
	  delete person;
	  person = sample[isample->second];
	  old_person++;
	}
      else
	new_person++;



      // Only look at phenotype if not already created

      if (!already_in)
	{
	  	
  	  //////////////////
	  // A non-founder?
	
	  person->founder = (person->pat == "0" && person->mat == "0") ? true : false;


	  /////////////////////////////////////////////////////
	  // Set missing status; test for quantitative traits?

	  if (phenotype == par::missing_phenotype)
	    person->missing = true;
	  else
	    {
	      if ( ! from_string<double>( person->phenotype, phenotype, std::dec ) )
		person->missing = true;
	      else
		if (phenotype != "0" && 
		    phenotype != "1" && 
		    phenotype != "2" ) 
		  {
		    par::qt = true;
		    par::bt = false; 
		  }
	    }
	  
	}
      

      /////////////////////////////////////// 
      // Add necessary space for a new person 
      // Missing genotypes by default

      if (!already_in)      
	{
	  if (par::SNP_major)
	    {
	      // Add a new missing person to each SNP
	      vector<CSNP*>::iterator s = SNP.begin();
	      while ( s != SNP.end() ) 
		{
		  (*s)->one.push_back(true);
		  (*s)->two.push_back(false);
		  s++;
		}	      
	      // And set the individual number
	      indn = n + new_person - 1;
	    }
	  else
	    {
	      // Add all new SNPs to this person
	      person->one.resize(locus.size(),true);
	      person->two.resize(locus.size(),false);
	    }

	}


      // Record whether this individual is new or not
      existing_person_list.push_back(already_in);	      
      sample2.push_back(person);
      
      // Add individual to list, if need be
      if (!existing_person_list[c])
	{
	  sample.push_back(person);
	  msample.insert(make_pair( person->fid+"_"+person->iid, 
				    sample.size()-1 ) );
	}
      
      // Increase person counter
      c++;
      
      
      
   }
  
  if (!par::merge_list)
    {
      printLOG(int2str( existing_person_list.size() ) 
	       + " individuals merged from [ " 
	       + par::merge_famfile + " ] \n");
      printLOG("Of these, " + int2str(new_person) 
	       + " were new, " + int2str( old_person ) +
	       " were already in current data\n\n");
    }




  ////////////////////////////////////
  // Read genotype information, merge

  ifstream BIT;
  bool bfile_SNP_major = openBinaryFile( par::merge_bedfile, BIT ); 
  
  if (bfile_SNP_major != par::SNP_major)
    error("BED files must both be SNP-major or both individual-major for merging\n");

  if ( (!par::SNP_major) || (!bfile_SNP_major) ) 
    error("Cannot --bmerge individual-major BED files -- convert to SNP-major");
  
  
   ///////////////////////////
   // SNP-major mode
  
   if (bfile_SNP_major)
     {
       
       // Person look-up table
       vector<Individual*>::iterator person = sample2.begin();
       vector<int> vindn;
       while ( person != sample2.end() )
	 {
	   map<string,int>::iterator isample = 
	     msample.find((*person)->fid+"_"+(*person)->iid);
	   int indn;
	   if ( isample != msample.end() ) 
	     indn = isample->second;
	   else
	     error("Internal error in --bmerge... should not happen...\n");	   
	   vindn.push_back(indn);
	   ++person;
	 }
       

      CSNP * snp;
      
      // Outer loop for SNPs
      int s=0;

      while (s<locus2.size()) // for all SNPs in second file
	{
	  
	  // Look up SNP information
	  
// 	  map<string,int>::iterator ilocus = mlocus.find(locus2[ s ]->name);
// 	  int l2 = ilocus->second;
	  
	  int k0 = (int)ordered[s].freq;
	  ilocus = mlocus.find(locus2[k0]->name);
	  int k = ilocus->second;
	  bool flipmode = flip_alleles.find( locus2[k0]->name ) != flip_alleles.end();
	  string snp1 = locus2[ k0 ]->allele1;
	  string snp2 = locus2[ k0 ]->allele2;
	  bool existence = exists.find(locus2[k0]->name)->second;
	  
	  
	  //////////////////////////////////////
	  // Inner loop for individuals
	  
	  vector<int>::iterator indn = vindn.begin();
	  
	  while ( indn != vindn.end() )
	    {
	      
	      char ch[1];
	      BIT.read(ch,1);
	      
	      bitset<8> b;
	      b = ch[0];	  

	      int c=0;
	      
	      if (!BIT) 
		error("Problem with the BED file...has the FAM file been changed?\n");
	      	      
	      while (c<7 && indn != vindn.end() ) 
		{
		  
		  bool s1 = b[c++];
		  bool s2 = b[c++];
		  
		  if ( flipmode && s1 == s2 )
		    {
		      s1 = !s1;
		      s2 = !s2;
		    }
		  
		  string one = snp1;
		  string two = snp2;
		  if (s1 && s2) one=snp2;
		  else if ( (!s1) && (!s2) ) two=snp1;
		  else if ( s1 && !s2 ) one=two=par::missing_genotype;
		  
		  bool already_in = false;
		  if ( *indn < original_sample_size ) already_in = true;
		  else already_in = existing_person_list[*indn - original_sample_size];

                  bool e = reconcileMerge( *indn, k , 
 					   one, two, 
	 				   already_in,
					   existence,
					   MERD, misstrand_dummy);	
		  
		  if (e) fatal_error=true;

		  // next person
		  indn++;		  
		  
		}
	      
	    }
	  
	  
	  // next SNP
	  s++;
	}
      
      // Set file mode
      par::SNP_major = true;
      
    }
   
   
   ////////////////////////////////////
   // Individual-major mode

//    else
//      {
       
//        // Outer loop for individuals
//        vector<Individual*>::iterator person = sample.begin();
//        while ( person != sample.end() )
// 	 {
	   
// 	   // Inner loop for SNPs
// 	   int s=0;
// 	   while (s<locus.size()) // for all SNPs
// 	     {
	       
// 	       char ch[1];
// 	       BIT.read(ch,1);
// 	       if (!BIT) error("Problem with the BED file... has the FAM file been changed?\n");
	       
// 	       bitset<8> b;
// 	       b = ch[0];	  
	       
// 	       int c=0;
	       
// 	       while (c<7 && s<locus.size())
// 		 {
// 		   (*person)->one[ s ] = b[c++];
// 		   (*person)->two[ s ] = b[c++];	      
// 		   s++;
// 		 }	 	  
// 	     }
	   
// 	   person++;
// 	 }
       
//        // Set file mode
//        par::SNP_major = false;
//      }
   
   
   
   
   
   //////////////
   // Close files
   
   FAM.clear();
   FAM.close();
   BIT.clear();
   BIT.close();
   
   
  
  // If a binary trait, now make 0 missing also
  // i.e. if we never saw other than missing, 0, 1 or 2
  
  if (par::bt)
    for (int i=0; i<sample.size(); i++)
      if ( sample[i]->phenotype == 0 )
	sample[i]->missing = true;


  
  if (par::merge_mode >=6) 
    {
      printLOG("Results from diff ( merge mode " 
	       + int2str(par::merge_mode) + " ) written to [ " +
	       par::output_file_name + ".diff ]\n");
      MERD.close();
      
      printLOG("Of " + int2str(diff_overlap) + " overlapping SNPs, " 
	       + int2str( diff_nonmissing_overlap ) + " were both genotyped\nand "
	       + int2str( diff_concordant_overlap ) + " were concordant\n");
      printLOG("Concordance rate is " 
	       + dbl2str( (double)diff_concordant_overlap 
			  / (double)diff_nonmissing_overlap ) + "\n");
      


      shutdown();
    }
  
  


  // Phenotype statistics
  if (!par::merge_list)
    {
      int nm=0;
      for (int i=0;i<sample.size();i++)
	if(!sample[i]->missing) nm++;
      printLOG(int2str( nm ) + " individuals with nonmissing phenotypes\n");
      
      if (par::bt) 
	{
	  printLOG("Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)\n");
	  if (par::missing_phenotype!="0")
	    printLOG("Missing phenotype value is also " + par::missing_phenotype + "\n");
	  
	  int ncase = 0;
	  int ncontrol = 0;
	  for (int i=0; i<sample.size(); i++)
	    if ( sample[i]->phenotype == 1 )
	      ncontrol++;
	    else if ( sample[i]->phenotype == 2 )
	      ncase++;
	  printLOG(int2str(ncase)+" cases and "+int2str(ncontrol)+" controls\n");
	}
      else 
	{
	  printLOG("Assuming a quantitative trait\n");
	  printLOG("Missing phenotype value is " + par::missing_phenotype + "\n");
	}
      
    }
  
}





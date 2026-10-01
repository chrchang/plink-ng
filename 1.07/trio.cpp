

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
#include <map>
#include <vector>
#include <set>
#include <cmath>

#include "plink.h"
#include "options.h"
#include "helper.h"
#include "crandom.h"
#include "sets.h"
#include "perm.h"

extern Plink * PP;

string gprint(int l, bool s1, bool s2);

////////////////////////////////////////////////////////
// Helper function: add individual to family, w/ checks.

void addParent(Family * f, Individual * person)
{
  // A real person? 
  if ( person )
    {
      // as father
      if ( person->sex )
	{
	  if ( f->pat == NULL ) { 
	    f->pat = person; 

	    // Do not overwrite (*see note below)
	    if (!person->family)
	      person->family = f; 
	  } 
	  else error("Problem pedigree structure: two fathers found : family "
		     +person->fid);
	  
	}
      else if (person->sexcode=="2")
	{   
	  // as mother
	  if ( f->mat == NULL ) { 
	    f->mat = person; 

	    // Do not overwrite (*see note below)
	    if (!person->family)
	      person->family = f; 
	  }
	  else error("Problem pedigree structure: two mothers found: family "
		     +person->fid);
	}
      else
	{
	  error("Problem with ambiguous parental sex codes for family "
		+person->fid);
	}
    }
  else
    {
      // Otherwise, add as a dummy individual
      if ( !f->pat )
	f->pat = person;
      else if ( ! f->mat )
	f->mat = person;
      else error("Internal error: allocated too many parents...\n");
    }
      
}

void addPerson(Family * f, Individual * person)
{
  
  // Add as child if does not already exist
  
  for (int c=0; c<f->kid.size(); c++)
    if ( f->kid[c]->iid == person->iid ) 
      error("Problem with family "
	    +f->kid[c]->fid+" child "+f->kid[c]->iid
	    +": offspring already exists\n");
  f->kid.push_back(person);
  // Always set (*see note)
  person->family = f;
}


// Note on priority of setting person->family pointer:
// This will be preferentially set to the family in which 
// the individual is the offspring (i.e. an individual can 
// only appear in one nuclear family as an offspring, but 
// multiple as a founder). Therefore, if nothing has been set 
// (singleton) set the family pointer for founder (i.e. so as 
// not to overwrite offspring family pointer) but otherwise
// always overwrite if setting the offspring family.

// In qfam.cpp then all family[] are used to construct B scores
// (i.e. an individual's genotypes might be used in several families)
// but the person->family is used to enter the person into the actual
// analysis (i.e. so the individual will always appear as offspring,
// or sibship (even is S=1)


void Plink::parseTrios()
{  

  /////////////////////////////////////
  // General check for unique FID, IIDs
  // and that no IID == 0 

  set<string> fid_iid;
  vector<Individual*>::iterator person = sample.begin();
  while ( person != sample.end() )
    {
      if ( (*person)->iid == "0" )
	error("Family "+(*person)->fid+" has person with reserved 0 ID\n");
      
      if ( (*person)->iid == (*person)->pat )
	error("Family "+(*person)->fid+" has person "+(*person)->iid+" who is own father");

      if ( (*person)->iid == (*person)->mat )
	error("Family "+(*person)->fid+" has person "+(*person)->iid+" who is own mother");

      string s = (*person)->fid+"_"+(*person)->iid;
      if ( fid_iid.find(s) != fid_iid.end() )
	error("Duplicate individual in pedigree: "
	      +(*person)->fid+" "
	      +(*person)->iid+"\n");
      else
	fid_iid.insert(s);
      person++;
    }
  


  /////////////////////////////////////
  // First consider all nonfounders 
  // with 2 parents

  map<string,Family*> fam;
  map<string,Family*>::iterator f;
  set<Individual*> infamily;
  
  person = sample.begin();
  
  while ( person != sample.end() )
    {
      
      // For non-founders
      if ( ! (*person)->founder ) 
	{
	  
	  string fpm = (*person)->fid+
	    "_"+(*person)->pat+
	    "_"+(*person)->mat;
	  
	  f = fam.find(fpm);
	  
	  // Have we already come across this parental pairing?
	  if ( f != fam.end() ) 
	    {
	      addPerson(f->second,*person);
	      infamily.insert(*person);
	    }
	  else 
	    {
	      Family * nfam = new Family;
	      
	      // Add person to new family
	      addPerson(nfam,*person);

	      infamily.insert(*person);

	      // And the parental pairing
	      
	      addParent(nfam,(*person)->pp);
	      addParent(nfam,(*person)->pm);		  
	      
	      infamily.insert( (*person)->pp );
	      infamily.insert( (*person)->pm );

	      if ( (*person)->pp != NULL && 
		   (*person)->pm != NULL )
		nfam->parents = true;
	      else
		nfam->sibship = true;
	      
	      // And add it to the list
	      fam.insert(make_pair(fpm,nfam));	  
	      
	    }
	}
      
      
      // Consider next individual
      person++;
    }
  
  // Rescan list for singleton founders
  person = sample.begin();
  while ( person != sample.end() )
    {
      if ( infamily.find( *person ) == infamily.end() )
	{
	  Family * nfam = new Family;
	  
	  // Have we already seen a singleton founder in this family?
	  map<string,Family*>::iterator myf = fam.find( (*person)->fid+"_0_0" );
	  if ( myf != fam.end() )
	    nfam = myf->second;
	      
	  addPerson(nfam,*person);
	  nfam->singleton = true;
	  fam.insert(make_pair((*person)->fid+"_0_0",nfam));	  
	}
      person++;
    }  
  

  // Counts
  int disc_parent_cnt = 0;
  int with_parents_ind_cnt = 0;
  int with_parents_fam_cnt = 0;
  int aff_with_parents_trio_cnt = 0;
  int without_parents_ind_cnt = 0;
  int without_parents_fam_cnt = 0;
  int singleton_cnt = 0;

  // Assign family type flags
  
  for ( f = fam.begin() ; f != fam.end() ; f++)
    {
      
      Family * mf = f->second;
      
      // TDT = 2 parents + atleast 1 affected offspring
      
      if ( mf->kid.size()>=0) 
	{
	  
	  if ( mf->parents )
	    {
	      with_parents_fam_cnt++;
	      for (int k=0; k<mf->kid.size(); k++)
		{
		  with_parents_ind_cnt++;
		  if (mf->kid[k]->aff)
		    {
		      aff_with_parents_trio_cnt++;
		      mf->TDT = true;
		    }
		}
	      
	    }
	  else
	    {
	      without_parents_fam_cnt++;
	      for (int k=0; k<mf->kid.size(); k++)
		without_parents_ind_cnt++;
	    }
	  
	  // Set flag for phenotypically discordant parents
	  
	  mf->discordant_parents = false;
	  if ( mf->parents )
	    if ( mf->pat->phenotype != mf->mat->phenotype &&
		 (!mf->pat->missing)  && 
		 (!mf->mat->missing)  )
	      {
		mf->discordant_parents = true;
		disc_parent_cnt++;
	      }
	} 
      
      if (mf->singleton) 
	singleton_cnt++;

      family.push_back(mf);
      
    }
  
  

  


  
  // Report counts
  
  printLOG(int2str(family.size())+" nuclear families, ");
  printLOG(int2str(singleton_cnt)+" founder singletons found\n");
  printLOG(int2str(with_parents_ind_cnt)+" non-founders with 2 parents in "
	   +int2str(with_parents_fam_cnt)+" nuclear families\n");
  printLOG(int2str(without_parents_ind_cnt-singleton_cnt)+" non-founders without 2 parents in "
	   +int2str(without_parents_fam_cnt-singleton_cnt)+" nuclear families\n");

  if (par::bt)
    printLOG(int2str(aff_with_parents_trio_cnt)+" affected offspring trios\n");
	      
  printLOG(int2str(disc_parent_cnt)+" phenotypically discordant parent pairs found\n");
  
  if (disc_parent_cnt>0) par::discordant_parents = true;

  // DFAM routine has it's own dump pedigree code

  if ( par::dumpped && ! par::sibTDT_test ) 
    {
      string str = par::output_file_name + ".pdump";
      printLOG("Dumping pedigree information to [ " + str + " ]\n");
      ofstream PD(str.c_str(),ios::out);

      vector<Family*>::iterator f = family.begin();
     
      while ( f != family.end() )
	{

	  if ( (*f)->singleton )
	    {
	      PD << "SINGLETON(S)\t" 
		 << (*f)->kid[0]->fid << " : ";
		for (int k=0; k < (*f)->kid.size() ;k++)
		  PD << (*f)->kid[k]->iid << " ";
	      PD << "\n";
	    }
	  else if ( (*f)->sibship )
	    {
	      PD << "SIBSHIP  \t" << (*f)->kid[0]->fid << " : ";
	      for ( int k=0; k<(*f)->kid.size(); k++)
		PD << (*f)->kid[k]->iid << " ";
	      PD << "\n";
	    }
	  else if ( (*f)->parents )
	    {
	      PD << "W/PARENTS\t" << (*f)->pat->fid << " : ";
	      PD << (*f)->pat->iid << " x " << (*f)->mat->iid << " -> ";
	      for ( int k=0; k<(*f)->kid.size(); k++)
		PD << (*f)->kid[k]->iid << " ";
	      PD << "\n";	     
	    }
	  else
	    PD << "UNDEFINED\t" 
		 << (*f)->pat->fid << " "
		 << (*f)->pat->iid << "\n";
	  
	  // Next family
	  f++;
	}
      
      PD << "\n\n";
      
      PD << "Listing by individual: columns are (0/1 for true false) \n"
	 << "  FID\tFamily ID\n"
	 << "  IID\tIndividual ID\n"
	 << "  Phenotype\n"
	 << "  Parents?\n"
	 << "  Singleton?\n"
	 << "  Sibship?\n"
	 << "  Discordant parents?\n"
	 << "  TDT?\n\n";
      
      // Now list by individual 
      
      vector<Individual*>::iterator person = sample.begin();

      while ( person != sample.end() ) 
	{
	  
	  PD << (*person)->fid << "\t"
	     << (*person)->iid << "\t"
	     << (*person)->phenotype << "\t"
	     << (*person)->family->parents << " "
	     << (*person)->family->singleton << " "
	     << (*person)->family->sibship << " " 
	     << (*person)->family->discordant_parents << " " 
	     << (*person)->family->TDT << "\n";
	    	    
	  person++;
	}

      PD.close();
    }
  

}



void Plink::checkMendel()
{


  //////////////////////////////////
  // Individual-major mode analysis
  
  if (par::SNP_major) SNP2Ind();

  ofstream MEN;
  ofstream MENL;
  ofstream MENI;
  ofstream MENF;

  if (par::MENDEL_test)
    printLOG("Filtering SNPs/families for Mendel Error rates above "+
	     dbl2str(par::MENDEL_snp)+", "+dbl2str(par::MENDEL_ind)+"\n");
  

  if (par::MENDEL_report)
    {
      string f = par::output_file_name+".mendel";
      string fi = par::output_file_name+".imendel";
      string ff = par::output_file_name+".fmendel";
      string fl = par::output_file_name+".lmendel";
      
      printLOG("Writing all Mendel errors to [ " + f +" ]\n");
      printLOG("Writing per-offspring Mendel summary to [ " + fi + " ]\n");
      printLOG("Writing per-family Mendel summary to [ " + ff + " ]\n");
      printLOG("Writing per-locus Mendel summary to [ " + fl + " ]\n");
      
      MEN.open(f.c_str(),ios::out);
      MENL.open(fl.c_str(),ios::out);
      MENI.open(fi.c_str(),ios::out);
      MENF.open(ff.c_str(),ios::out);
      
      MEN << setw(par::pp_maxfid) << "FID" << " " 
	  << setw(par::pp_maxiid) << "KID" << " " 
	  << setw(4) << "CHR" << " " 
	  << setw(par::pp_maxsnp) << "SNP" << " " 
	  << setw(6) << "CODE"   
	  << setw(22) << "ERROR"  
	  << "\n";
      
      MENL << setw(4) << "CHR" << " " 
	   << setw(par::pp_maxsnp) << "SNP" << " " 
	   << setw(4) << "N\n";
      
      MENI << setw(par::pp_maxfid) << "FID"  << " " 
	   << setw(par::pp_maxiid) << "IID" << " " 
	   << setw(4) << "N\n";
      
      MENF << setw(par::pp_maxfid) << "FID" << " " 
	   << setw(par::pp_maxiid) << "PAT" << " "
	   << setw(par::pp_maxiid) << "MAT" << " "
	   << setw(6) << "CHLD" << " " 
	   << setw(4) << "N" << "\n";
    }


  // Flag to remove SNP if need be
  vector<bool> mendel_locus(nl_all,false);
  // Counts per family
  vector<int> mendel_family(family.size(),0);
  vector<int> mendel_pat(family.size(),0);
  vector<int> mendel_mat(family.size(),0);
  vector<vector<int> > mendel_indiv(family.size());
  for (int f=0; f<family.size(); f++) 
    mendel_indiv[f].resize(family[f]->kid.size());
  
  // Count number of trios
  int n_trios = 0; 
  for (int f=0; f<family.size(); f++) 
    n_trios += family[f]->kid.size();
  
  int total_m = 0; // total number of Mendel errors
  int l_removed = 0; // # of SNPs removed due to high ME rate

  // Test each locus
  for (int l=0; l<nl_all; l++)
    {
      
      // Skip haploid markers
      if (par::chr_haploid[locus[l]->chr]) continue;
      
      // Flag for X markers
      bool X = par::chr_sex[locus[l]->chr];
      
      int m=0; // number of Mendel errors per locus
      
      for (int f=0; f<family.size(); f++) 
	{
	  
	  // Have we observed parents?
	  if ( ! family[f]->parents ) continue;
	  
	  Individual * pat = family[f]->pat;
	  Individual * mat = family[f]->mat;
	  vector<Individual *> kid = family[f]->kid;
	  
	  mendel_indiv[f].resize(kid.size());

	  for (int c=0; c<kid.size(); c++)
	    {
	      
	      // If child genotype is missing, skip
	      if ( kid[c]->one[l] && !kid[c]->two[l] ) continue;

	      int mendel_type = 0;
	      
	      // For autosomal markers
	      if ( (!X) || (!kid[c]->sex) )
		{
		  
		  if ( (!kid[c]->one[l]) && kid[c]->two[l] )
		    {
		      
		      // KID = 01
		      // 00x00 -> 01  (m1)
		      // 11x11 -> 01  (m2)
		      
		      
		      if (  ( (!pat->one[l]) && (!pat->two[l]) ) && 
			    ( (!mat->one[l]) && (!mat->two[l]) ) )
			{
			  if ( ! par::preserve_mendel_errors)
			    {
			      kid[c]->one[l] = true;
			      kid[c]->two[l] = false;

			      mat->one[l] = true;
			      mat->two[l] = false;

			      pat->one[l] = true;
			      pat->two[l] = false;
			    }

			  m++;
			  mendel_type = 1;
			}
		      else if (  ( pat->one[l] && pat->two[l] ) && 
				 ( mat->one[l] && mat->two[l] ) )
			{
			  if ( ! par::preserve_mendel_errors)
			    {
			      kid[c]->one[l] = true;
			      kid[c]->two[l] = false;
			      
			      mat->one[l] = true;
			      mat->two[l] = false;
			      
			      pat->one[l] = true;
			      pat->two[l] = false;
			    }
			  m++;
			  mendel_type = 2;
			}
		    }
		  else if ( (!kid[c]->one[l]) && (!kid[c]->two[l]) )
		    {
		      // KID = 00
		      // 00x11 -> 00 (m3) P11->00
		      // 01x11 -> 00 (m3) 
		      // ??x11 -> 00 (m3)
		      
		      // 11x00 -> 00 (m4) M11->00
		      // 11x01 -> 00 (m4) 
		      // 11x?? -> 00 (m4) 
		      
		      // 11x11 -> 00 (m5) P11+M11->00 
		      
		      // Hom parent --> opposite hom child
		      
		      // rule = at least one '11' parent
		      
		      if ( ( pat->one[l] && pat->two[l] ) || 
			   ( mat->one[l] && mat->two[l] ) )
			{
			  if ( ! par::preserve_mendel_errors )
			    {
			      kid[c]->one[l] = true;
			      kid[c]->two[l] = false;
			    }
			  m++;
			  
			  if ( pat->one[l] && pat->two[l] && 
			       mat->one[l] && mat->two[l] )
			    mendel_type = 5;
			  else if ( pat->one[l] && pat->two[l] )
			    {
			      mendel_type = 3;
			      
			      if ( ! par::preserve_mendel_errors )
				{
				  pat->one[l] = true;
				  pat->two[l] = false;			      
				}
			    }
			  else
			    {
			      mendel_type = 4;		      
			      if ( ! par::preserve_mendel_errors )
				{
				  mat->one[l] = true;
				  mat->two[l] = false;
				}
			    }
			  
			}
		    }	
		  else 
		    {
		      
		      // KID = 11
		      
		      // 00x01 -> 11 (m6)
		      // 00x11 -> 11
		      // 00x?? -> 11
		      
		      // 01x00 -> 11 (m7)
		      // 11x00 -> 11
		      // ??x00 -> 11
		      
		      // 00x00 -> 11 (m8) P00+M00->11		 
		      
		      // rule = at least one '00' parent
		      
		      if ( ( (!pat->one[l]) && (!pat->two[l]) ) || 
			   ( (!mat->one[l]) && (!mat->two[l]) ) )
			{
			  if ( ! par::preserve_mendel_errors)
			    {
			      kid[c]->one[l] = true;
			      kid[c]->two[l] = false;
			    }
			  m++;
			  
			  if ( (!pat->one[l]) && (!pat->two[l]) && 
			       (!mat->one[l]) && (!mat->two[l]) )
			    mendel_type = 8;
			  else if ( (!pat->one[l]) && (!pat->two[l]) )
			    {
			      mendel_type = 6;
			      if ( ! par::preserve_mendel_errors )
				{
				  pat->one[l] = true;
				  pat->two[l] = false;			     
				}
			    }
			  else
			    {
			      mendel_type = 7;		      
			      if ( ! par::preserve_mendel_errors )
				{
				  mat->one[l] = true;
				  mat->two[l] = false;
				}
			    }
			  
			}
		      		      
		    }
		  
		}
	      else
		{
		  if ( kid[c]->one[l] && kid[c]->two[l] && 
		       (!mat->one[l]) && (!mat->two[l]) )
		    {
		      m++;
		      mendel_type = 9;
		      if ( ! par::preserve_mendel_errors)
			{
			  kid[c]->one[l] = true;
			  kid[c]->two[l] = false;
			  
			  mat->one[l] = true;
			  mat->two[l] = false;
			}
		    }
		  
		  if ( (!kid[c]->one[l]) && (!kid[c]->two[l]) && 
		       mat->one[l] && mat->two[l] )
		    {
		      m++;
		      mendel_type = 10;
		      if ( ! par::preserve_mendel_errors)
			{
			  kid[c]->one[l] = true;
			  kid[c]->two[l] = false;		      
			  
			  mat->one[l] = true;
			  mat->two[l] = false;
			}
		    }
		  
		}
		    
	      // Individual counts
	      
	      // m1   00x00  ->  01   K / P+M
	      // m2   11x11  ->  01   K / P+M

	      // m3   11x**  ->  00   K / P
	      // m4   **x11  ->  00   K / M
	      // m5   11x11  ->  00   K

	      // m6   00x**  ->  11   K / P
	      // m7   **x00  ->  11   K / M
	      // m8   00x00  ->  11   K

	      // X marker errors for male offspring
	      // m9   **x00 -> 11     K / M
	      // m10  **x11 -> 00     K / M
	      

	      if (par::MENDEL_report && mendel_type>0)
		{
		  MEN << setw(par::pp_maxfid) << family[f]->pat->fid  << " " 
		      << setw(par::pp_maxiid) << kid[c]->iid << " " 
		      << setw(4)  << locus[l]->chr  << " " 
		      << setw(par::pp_maxsnp) << locus[l]->name  << " " ;
		  
		  if (mendel_type==1) 
		    {
		      MEN << setw(6) << "1"<< setw(22);
		      string s = locus[l]->allele1 + "/" + locus[l]->allele1 + " x " 
			+ locus[l]->allele1 + "/" + locus[l]->allele1 + " -> " 
			+ locus[l]->allele1 + "/" + locus[l]->allele2;
		      MEN << s << "\n";
		    }
		  else if (mendel_type==2)
		    {
		      MEN << setw(6) << "2"<< setw(22);
		      string s = locus[l]->allele2 + "/" + locus[l]->allele2 + " x "
			+ locus[l]->allele2 + "/" + locus[l]->allele2 + " -> "
			+ locus[l]->allele1 + "/" + locus[l]->allele2;
		      MEN << s << "\n";
		    }
		  else if (mendel_type==3)
		    {
		      MEN << setw(6) << "3"<< setw(22);
		      string s = locus[l]->allele2 + "/" + locus[l]->allele2 + " x "
			+ "*/*" + " -> "
			+ locus[l]->allele1 + "/" + locus[l]->allele1;
		      MEN << s << "\n";
		    } 
		  else if (mendel_type==4)
		    {	
		      MEN << setw(6) << "4"<< setw(22);
		      string s = string("*/*") + string(" x ") 
			+ locus[l]->allele2 + "/" + locus[l]->allele2 + " -> "
			+ locus[l]->allele1 + "/" + locus[l]->allele1;
		      MEN << s << "\n";
		    }	
		  else if (mendel_type==5)
			{	
			  MEN << setw(6) << "5"<< setw(22);
			  string s= locus[l]->allele2 + "/" + locus[l]->allele2 + " x "
			    + locus[l]->allele2 + "/" + locus[l]->allele2 + " -> "
			    + locus[l]->allele1 + "/" + locus[l]->allele1;
			  MEN << s << "\n";
			}	
		  else if (mendel_type==6)
		    {
		      MEN << setw(6) << "6"<< setw(22);
		      string s = locus[l]->allele1 + "/" + locus[l]->allele1 + " x "
 			+ "*/*" + " -> "
			+ locus[l]->allele2 + "/" + locus[l]->allele2;
		      MEN << s << "\n";
		    }
		  else if (mendel_type==7)
		    {
		      MEN << setw(6) << "7"<< setw(22);
		      string s = string("*/*") + string(" x ") 
			+ locus[l]->allele1 + "/" + locus[l]->allele1 + " -> "
			+ locus[l]->allele2 + "/" + locus[l]->allele2;
			MEN << s << "\n";
		    }
		  else if (mendel_type==8)
		    {
		      MEN << setw(6) << "8" << setw(22); 
		      string s = locus[l]->allele1 + "/" + locus[l]->allele1 + " x "
			+ locus[l]->allele1 + "/" + locus[l]->allele1 + " -> "
			+ locus[l]->allele2 + "/" + locus[l]->allele2;
			MEN << s << "\n";
		    }
		  else if (mendel_type==9)
		    {
		      MEN << setw(6) << "9" << setw(22); 
		      string s = string("*/*") + " x "
			+ locus[l]->allele1 + "/" + locus[l]->allele1 + " -> "
			+ locus[l]->allele2 + "/" + locus[l]->allele2;
		      MEN << s << "\n";
		    }
		  else if (mendel_type==10)
		    {
		      MEN << setw(6) << "10" << setw(22); 
		      string s = string("*/*") + " x "
			+ locus[l]->allele2 + "/" + locus[l]->allele2 + " -> "
			+ locus[l]->allele1 + "/" + locus[l]->allele1;
		      MEN << s << "\n";
		    }

		}
	      
	      
	      // Family count
	      if (mendel_type>0) 
		{
		  // We may wish to add other weighting schemes here
		  
		  if ( mendel_type==1 || mendel_type==2 ) 
		    { 
		      mendel_indiv[f][c]++; 
		      mendel_pat[f]++; 
		      mendel_mat[f]++; 
		    }
		  else if ( mendel_type==5 || mendel_type==8 ) 
		    mendel_indiv[f][c]++;
		  else if ( mendel_type==3 || mendel_type==6 ) 
		    { 
		      mendel_indiv[f][c]++; 
		      mendel_pat[f]++; 
		    }
		  else if ( mendel_type==4 || mendel_type==7 ) 
		    { 
		      mendel_indiv[f][c]++; 
		      mendel_mat[f]++; 
		    }	      
		  else if ( mendel_type==9 || mendel_type==10 ) 
		    { 
		      mendel_indiv[f][c]++; 
		      mendel_mat[f]++; 
		    }
		}

	      // Family count
	      if (mendel_type>0) 
		mendel_family[f]++;
	      
	    }

	}


      // Report per-SNP Mendel count 
      if (par::MENDEL_report)
	MENL << setw(4)  << locus[l]->chr << " "
	     << setw(par::pp_maxsnp) << locus[l]->name << " "
	     << setw(4)  << m << "\n";

      // Is this a bad SNP? 
      if ( (double)m/(double)n_trios > par::MENDEL_snp)
	{
	  l_removed++;
	  mendel_locus[l] = true;
	}

      // Keep track of total number of errors
      total_m += m;
      
    }


  printLOG(int2str(total_m)+" Mendel errors detected in total\n");

  if (par::MENDEL_report)
    {  
      for (int f=0; f<family.size(); f++) 
	{
	  if ( family[f]->parents )
	    MENF << setw(par::pp_maxfid) << family[f]->pat->fid << " "
		 << setw(par::pp_maxiid) << family[f]->pat->iid << " "
		 << setw(par::pp_maxiid) << family[f]->mat->iid << " "
		 << setw(6) << family[f]->kid.size() << " "
		 << setw(4) << mendel_family[f]  
		 << "\n";
	}
      
      for (int f=0; f<family.size(); f++) 
	{
	  if ( family[f]->parents )
	    {
	      vector<Individual *> kid = family[f]->kid;
	      MENI << setw(par::pp_maxfid) << family[f]->pat->fid << " "
		   << setw(par::pp_maxiid) << family[f]->pat->iid << " "
		   << setw(4)  << mendel_pat[f] 
		   << "\n";
	      MENI << setw(par::pp_maxfid) << family[f]->mat->fid << " "
		   << setw(par::pp_maxiid) << family[f]->mat->iid << " "
		   << setw(4)  << mendel_mat[f]  << " "
		   << "\n";
	      for (int c=0; c<kid.size(); c++)
		MENI << setw(par::pp_maxfid) << family[f]->pat->fid << " "
		     << setw(par::pp_maxiid) << kid[c]->iid << " "
		     << setw(4)  << mendel_indiv[f][c] 
		     << "\n";
	    }      
	}

      MEN.close();
      MENL.close();
      MENI.close();
      MENF.close();

      shutdown();
    }



  // Using Mendel error rates to automatically remove SNPs 
  // and families
  
  if (par::MENDEL_test)
    {
      
      ////////////////////////////////////////
      // Remove selected loci from locus list, 
      // by copying rest to a new list
      // People/genotypes first, then locus/map info

      bool any_mendel = false;
      for (int l=0; l<nl_all; l++)
	if ( mendel_locus[l] ) any_mendel = true;
      
      int orig_nl = nl_all; // copy for comparison below

      if (any_mendel)
	int tmp = deleteSNPs(mendel_locus);
      

      map<Individual*,int> badfam;
      int f_removed = 0;
      for (int f=0;f<family.size();f++)
	{
	  if ((double)mendel_family[f]/(double)orig_nl > par::MENDEL_ind) 
	    {
	      f_removed++;
	      badfam.insert(make_pair(family[f]->pat,0));
	      badfam.insert(make_pair(family[f]->mat,0));
	      for (int c=0; c<family[f]->kid.size(); c++)
		badfam.insert(make_pair(family[f]->kid[c],0));
	    }
	}
      
      // Remove individuals as appropriate
      int n_removed = 0;
      vector<bool> indel(sample.size(),false);

      for (int i=0;i<n;i++)
	{
	  if ( badfam.find(sample[i]) != badfam.end() )
	    {
	      indel[i] = true;
	      n_removed++;
	    }
	}

      if (n_removed>0)
	int tmp = deleteIndividuals(indel);
      
      
      printLOG(int2str(f_removed)+" families ( "+
	       int2str(n_removed)+" individuals ) removed due to Mendel errors\n");
      printLOG(int2str(l_removed)+" markers removed due to Mendel errors, "+
	       int2str(nl_all)+" remaining\n");            


      ////////////////////////////
      // Rebuild family structure?
      
      if (f_removed+n_removed>0)
	{

	  // Wipe existing family structure

	  printLOG("Rebuilding families after filtering on Mendel errors\n");

	  family.clear();

	  map<string,Individual*> fnd;
	  map<Individual*,int> idmap;

	  linkRelateds(idmap, fnd);
	  parseTrios();      
	}

    }

  
}


void Plink::makeFounders()
{

  for (int i=0; i<n; i++)
    {
      Individual * person = sample[i];
      
      if ( ! person->founder )
	{
	  Individual * father = sample[i]->pp;
	  Individual * mother = sample[i]->pm;
	  
	  if ( ! ( father && mother ) )
	  {
	    person->pat = person->mat = "0";
	    ++cnt_f;
	  }
	       
	}
    }
}


void Plink::pseudoCaseControl()
{
  printLOG("Writing pseudo case/control units to [ "+par::output_file_name+".tucc.ped ]\n");
  
  // Consider each trio unit, 
  // then all SNPs
  
  vector<bool> t1(nl_all);
  vector<bool> t2(nl_all);

  vector<bool> u1(nl_all);
  vector<bool> u2(nl_all);

  if ( par::SNP_major ) 
    SNP2Ind();
  
  ofstream POUT;
  POUT.open( (par::output_file_name+".tucc.ped").c_str(), ios::out);
  
  for (int f=0; f<family.size(); f++) 
    {
      
      // Have we observed parents?
      if ( ! family[f]->parents ) continue;
      
      Individual * pat = family[f]->pat;
      Individual * mat = family[f]->mat;
      vector<Individual *> kid = family[f]->kid;

      for (int c=0; c<kid.size(); c++)
	{
	  
	  Individual * child = family[f]->kid[c];
	  
	  // Score for each SNP
	  
	  for (int l=0; l<nl_all; l++)
	    {
	      
	      // If any individual is missing, we skip entirely
	      bool mat1 = mat->one[l];
	      bool mat2 = mat->two[l];
	      
	      bool pat1 = pat->one[l];
	      bool pat2 = pat->two[l];
	      
	      bool kid1 = child->one[l];
	      bool kid2 = child->two[l];

	      
	      ////////////////////////
	      // Missing ?

	      if ( ( pat1 && ! pat2 ) ||
		   ( mat1 && ! mat2 ) ||
		   ( kid1 && ! kid2 ) ) 
		{
		  t1[l] = true;
		  t2[l] = false;		  
		  u1[l] = true;
		  u2[l] = false;
		}
	      else
		{
		  bool X = par::chr_haploid[locus[l]->chr];
		  bool haploid = par::chr_sex[locus[l]->chr];
		  
		  // Autosome
		  if ( X || haploid ) 
		    {
		      t1[l] = true;
		      t2[l] = false;		  
		      u1[l] = true;
		      u2[l] = false;
		    }
		  else
		    {
		      // Transmitted alleles
		      t1[l] = kid1;
		      t2[l] = kid2;
		      
		      // Untransmitted alleles
		      int aCount = 0;
		      if ( pat1 ) ++aCount;
		      if ( pat2 ) ++aCount;
		      if ( mat1 ) ++aCount;
		      if ( mat2 ) ++aCount;
		      if ( kid1 ) --aCount;
		      if ( kid2 ) --aCount;
		      
		      if ( aCount == 0 )
			{
			  u1[l] = false;
			  u2[l] = false;
			}
		      else if ( aCount == 1 )
			{
			  u1[l] = false;
			  u2[l] = true;
			}
		      else if ( aCount == 2 )
			{
			  u1[l] = true;
			  u2[l] = true;
			}
		      else
			{
			  cout << kid1 << kid2 << " <- "
			       << pat1 << pat2 << " " 
			       << mat2 << mat2 << "\n";

			  error("Internal problem in --tucc");
			}
		    }
		} // Next SNP
	    }
      
	  // Output two rows in PED file for this trio
	  POUT << child->fid << " " 
	       << child->iid << "_T 0 0 " 
	       << child->sexcode << " 2 ";
	  for (int l=0; l<nl_all; l++)
	    POUT << gprint(l,t1[l],t2[l]);
	  POUT << "\n";
	  
	  POUT << child->fid << " " 
	       << child->iid << "_U 0 0 " 
	       << child->sexcode << " 1 ";
	  for (int l=0; l<nl_all; l++)
	    POUT << gprint(l,u1[l],u2[l]);
	  POUT << "\n";
	  
	}
      
    }
  
  POUT.close();
  
}


string gprint(int l, bool s1, bool s2)
{
  string a1 = par::recode_12 ? "1" : PP->locus[l]->allele1;
  string a2 = par::recode_12 ? "2" : PP->locus[l]->allele2;
  if      ( (!s1) && (!s2) )  
    return par::recode_delimit+a1+par::recode_indelimit+a1;
  else if ( (!s1) && s2 ) 
    return par::recode_delimit+a1+par::recode_indelimit+a2;
  else if (  s1   && s2 ) 
    return par::recode_delimit+a2+par::recode_indelimit+a2;
  else 
    return par::recode_delimit + par::missing_genotype
      + par::recode_indelimit+par::missing_genotype;
  
  return "?";
}


void Plink::makeMissingParents()
{
  // Add to a separate pile of dummy people

  map<string,Individual*> padded;
  
  for (int i=0; i<n; i++)
    {
      Individual * person = sample[i];
      
      if ( person->pp == NULL && person->pat != "0" )
	{      

	  Individual * d; 
	  string pcode = person->fid+"_"+person->pat;
	  map<string,Individual*>::iterator f = padded.find( pcode );
	  
	  if ( f != padded.end() )
	    d = f->second;
	  else
	    {
	      d = new Individual;
	      d->fid = person->fid;
	      d->iid = person->pat;
	      d->sex = true;
	      padded.insert( make_pair( pcode, d ) );
	    }
	  person->pp = d;
	  d->kids.push_back(person);	  
	}

      if ( person->pm == NULL && person->mat != "0" )
	{      
	  Individual * d; 
	  string pcode = person->fid+"_"+person->mat;
	  map<string,Individual*>::iterator f = padded.find( pcode );
	  if ( f != padded.end() )
	    d = f->second;
	  else
	    {
	      d = new Individual;
	      d->fid = person->fid;
	      d->iid = person->mat;
	      d->sex = false;
	      padded.insert( make_pair( pcode, d ) );
	    }	  

	  person->pm = d;
	  d->kids.push_back(person);
	}
      
      // And make sure that the actual parents also list who is their
      // child (i.e. if one parent was previously missing, this might
      // not be the case)
      
      if (  person->pm )
	{
	  bool found = false;
	  for ( int k=0; k<person->pm->kids.size(); k++)
	    if ( person->pm->kids[k] == person )
	      found = true;
	  if ( ! found ) 
	    person->pm->kids.push_back( person );
	}
      if ( person->pp )
	{
	  bool found = false;
	  for ( int k=0; k<person->pp->kids.size(); k++)
	    if ( person->pp->kids[k] == person )
	      found = true;
	  if ( ! found ) 
	    person->pp->kids.push_back( person );
	}

    }
  
  if (padded.size() > 0) 
    printLOG("Added " + int2str(padded.size()) + " dummy parents\n");
}

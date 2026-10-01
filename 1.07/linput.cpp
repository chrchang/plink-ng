

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
#include <map>
#include <algorithm>
#include <bitset>
#include <limits>
#include <errno.h>

#include "plink.h"
#include "options.h"
#include "helper.h"

extern ofstream LOG;

void Plink::readDataLongFormat()
{

  
  //////////////////////
  // Check files exist
  
  checkFileExists(par::lpedfile);
  checkFileExists(par::mapfile);
  checkFileExists(par::famfile);


  ///////////////////////////////////////////////
  // .map file
  
  vector<bool> include;
  vector<int> include_pos(0);
  int nl_actual=0;
  
  // Read in MAP file: this function also allocates 
  // slots for SNP-major mode

  readMapFile(par::mapfile,
	      include,
	      include_pos,
	      nl_actual);
  
  

  //////////////////////////////////////////////
  // First read a reference file? 

  map<string,string> refallele;
  map<string,string> refallele2;

  if ( par::ref_file )
    {
      set<string> mset;
      for (int l=0; l< locus.size(); l++)
	mset.insert( locus[l]->name );

      int notfound = 0;
      checkFileExists( par::ref_file_name );
      ifstream REF( par::ref_file_name.c_str(), ios::in );
      while ( ! REF.eof() )
	{
	  vector<string> tok = tokenizeLine( REF );
	  if ( tok.size() == 0 ) 
	    continue;
	  
	  if ( tok.size() != 2 && tok.size() != 3 )
	    error("Problem with line in [ " 
		  + par::ref_file_name 
		  + " ] :\n " + displayLine( tok ) );


	  if ( mset.find( tok[0] ) == mset.end() )
	    {
	      ++notfound;
	      continue;
	    }
	  
	  refallele.insert( make_pair( tok[0], tok[1] ));

	  // A second allele also specified?
	  if ( tok.size() == 3 ) 
	    refallele2.insert( make_pair( tok[0], tok[2] ));

	}

      printLOG("Read reference alleles for " 
	       + int2str( refallele.size() ) 
	       + " sites\n");

      if ( notfound>0 )
	printLOG(int2str( notfound ) 
		 + " SNPs in reference but not in map file\n");
      if ( refallele.size() < locus.size() ) 
	printLOG(int2str( locus.size() - refallele.size() ) 
		 + " SNPs in map file but not in reference\n");
      REF.close();
    }


  ///////////////////////////////////////////////
  // .fam
  
  readFamFile(par::famfile);
  

  // Allocate space for individual-major mode, set to 
  // missing by default...

  // Either missing (TF) by default; or reference allele (FF)
  bool code = ! par::ref_file;

  if ( ! par::SNP_major)
    {
      for (int i=0; i<sample.size(); i++)
	{
	  sample[i]->one.resize(nl_actual,code);
	  sample[i]->two.resize(nl_actual,false);
	}
    }
  else
    {
      for (int l=0; l<SNP.size(); l++)
	{
	  SNP[l]->one.resize(sample.size(),code);
	  SNP[l]->two.resize(sample.size(),false);
	}
    }

  
  if ( par::ref_file )
    {
      for (int l=0; l< locus.size(); l++)
	{	  
	  map<string,string>::iterator i = refallele.find( locus[l]->name );
	  map<string,string>::iterator i2 = refallele2.find( locus[l]->name );
	  
	  // If we cannot find, we need to set genotypes to missing instead
	  
	  if ( i != refallele.end() )
	    {
	      locus[l]->allele1 = i->second;
	      
	      if ( i2 != refallele2.end() )
		locus[l]->allele2 = i2->second;
	      else if ( par::lfile_allele_count )
		locus[l]->allele2 = i->second + "v";
	    }
	  else
	    {
	      for (int i=0; i<sample.size(); i++)
		{
		  if ( par::SNP_major )
		    {
		      SNP[l]->one[i] = true;
		      SNP[l]->two[i] = false;
		    }
		  else
		    {
		      sample[i]->one[l] = true;
		      sample[i]->two[l] = false;
		    }
		}
	    }
	}
    }

  
  
  ///////////////////////////////////////////////
  // .lgen
  
  FILE * PED;
  
  PED = fopen64(par::lpedfile.c_str(),"r");
  if ( PED == NULL )
    error("Problem opening LGEN file, errno = "+int2str(errno));
  
  
  // We can now read any number of individual/genotype lines, in any
  // order; we also do not assume that all genotypes are given --
  // these will be missing by default

  map<string,int> imap;
  map<string,int> iperson;

  for (int i=0; i<include.size(); i++)
    {
      if ( include[i] ) 
	{
	  int k = include_pos[i];
	  imap.insert( make_pair ( locus[k]->name , k ) );
	}
    }
  
  for (int i=0; i<sample.size(); i++)
    {
      iperson.insert( make_pair 
		      ( sample[i]->fid + "_" + sample[i]->iid , 
			i ) );
    }
  
  
    // Whether or not we want to look at a locus is in the include[] vector
    // The genomic position of locus i is k=include_pos[i] -> locus[k]

    bool fatal = false;
    string fmsg = "";
  
    while( ! feof(PED) )
      {
	
	string fid = "";
	string iid = "";
	string snp = "";
	string one = "";
	string two = "";
	
	int f = 0;
	
	if ( readString( PED , fid ) ) f++;
	
	if ( fid == "" ) 
	  continue;

	if ( readString( PED , iid ) ) f++;
	if ( readString( PED , snp ) ) f++;
	if ( readString( PED , one ) ) f++;
	

	map<string,int>::iterator im = imap.find(snp);
	
	int k = 
	  im != imap.end() ?
	  im->second : -1;

	// Need to read second allele?

	if ( ! ( par::compound_genotype_code || par::lfile_allele_count ) )
	  {
	    if ( readString( PED , two ) ) f++;
	  }
	else
	  {
	    if ( par::compound_genotype_code )
	      {
		if ( one.size() != 2 )
		  error("Problem with compound genotype not of length 2: [ " + one + " ]");
		two = one[1];
		one = one[0];
	      }
	    else if ( par::lfile_allele_count && k != -1 )
	      {
		// expect either a 0,1 or 2,  or missing code (anything other than 0,1 or 2)
		int a;
		if ( ! from_string<int>( a, one, std::dec ) )
		  a = -1;
		if ( a < 0 || a > 2 ) 
		  {
		    one = two = par::missing_genotype;
		  }
		else if ( a == 1 )
		  {
		    one = locus[k]->allele1;
		    two = locus[k]->allele2;
		  }
		else if ( a == 2 )
		  one = two = locus[k]->allele2;
		else if ( a == 0 )
		  one = two = locus[k]->allele1;
	      }
	  }
      
//      	cout << f << " " << "[" << fid << "] "
//      	     << "[" << iid << "] "
//      	     << "[" << snp << "] "
//      	     << "[" << one << "] "
//      	     << "[" << two << "] \n";

	
	map<string,int>::iterator peri 
	  = iperson.find( fid+"_"+iid );      
	
	Individual * person = 
	  peri != iperson.end() ?
	  sample[peri->second] : NULL ; 
	
	
	// Ignore this genotype?

	if ( ( ! person ) || k < 0 ) 
	  continue;

	int ip = peri->second;	
	Locus * loc = locus[k];
	
	/////////////////////////////////////////
	// Add allele names to list, if needed
	
	// If allele is not missing...
	if (one!=par::missing_genotype && two!=par::missing_genotype)
	  {
	    // ...and not already listed
	    if (one!=loc->allele1 && one!=loc->allele2)
	      {
		// ...then add to first empty slot
		if(loc->allele1=="") loc->allele1=one;
		else if(loc->allele2=="") loc->allele2=one;
		else {
		  // .. or show an error if no empty slots
		  if (!fatal)
		    fmsg = "Locus " + 
		      loc->name + " has >2 alleles:\n       individual " 
		      + person->fid + " " + person->iid + " has genotype [ " + one +" "+two+" ]\n"
		      + "       but we've already seen [ " + loc->allele1 + " ] and [ " + loc->allele2 + " ]\n";
		  fatal=true;
		  
		}
	      }
	  }
	      
	// Repeat for second allele, if different
	if (two!=one)
	  {
	    // If allele is not missing...
	    if (one!=par::missing_genotype)
	      // ...and not already listed
	      if (two!=loc->allele1 && two!=loc->allele2)
		{
		  // ...then add to first empty slot
		  if(loc->allele1=="") loc->allele1=two;
		  else if(loc->allele2=="") loc->allele2=two;
		  else {
		    if (!fatal)
		      fmsg = "Locus " + 
			loc->name + " has >2 alleles:\n       individual " 
			+ person->fid + " " + person->iid + " has genotype [ " + one +" "+two+" ]\n"
			+ "       but we've already seen [ " + loc->allele1 + " ] and [ " + loc->allele2 + " ]\n";
		    fatal=true;
		  }
		}
	  }
	      

	// Give an error message

	if ( fatal )
	  error(fmsg);


	/////////////////////////////
	// Add specific genotypes
	
	if (par::SNP_major)
	  {
	    
	    // 00 hom
	    if (one==loc->allele1 && two==loc->allele1)
	      {
		SNP[k]->one[ip] = false;
		SNP[k]->two[ip] = false;
	      }
	    
	    
	    // 01 het
	    else if (one!=par::missing_genotype && 
		     two!=par::missing_genotype && 
		     one!=two)	    
	      {
		SNP[k]->one[ip] = false;
		SNP[k]->two[ip] = true;
	      }
	    
	    // 11 hom
	    else if (one==loc->allele2 && two==loc->allele2)
	      {
		SNP[k]->one[ip] = true;
		SNP[k]->two[ip] = true;
	      }
	    
	    // 10 missing
	    else if (one==par::missing_genotype || two==par::missing_genotype)
	      {
		SNP[k]->one[ip] = true;
		SNP[k]->two[ip] = false;
	      }
	  }
	else
	  {
	    
	    // 00 hom
	    if (one==loc->allele1 && two==loc->allele1)
	      {
		person->one[k]=false;
		person->two[k]=false;
	      }
	    
	    
	    // 01 het
	    else if (one!=par::missing_genotype && 
		     two!=par::missing_genotype && 
		     one!=two)	    
	      {
		person->one[k]=false;
		person->two[k]=true;		  
	      }
	    
	    // 11 hom
	    else if (one==loc->allele2 && two==loc->allele2)
	      {
		person->one[k]=true;
		person->two[k]=true;		  
	      }
	    
	    // 10 missing
	    else if (one==par::missing_genotype || two==par::missing_genotype)
	      {
		person->one[k]=true;
		person->two[k]=false;		  
	      }
	  }
      }
    

    
    fclose(PED);
  
}



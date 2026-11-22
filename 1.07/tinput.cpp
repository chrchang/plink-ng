

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
#include <limits>
#include <errno.h>
#include "plink.h"
#include "options.h"
#include "helper.h"

extern ofstream LOG;

void Plink::readTransposedData()
{
  
  //////////////////////
  // Check files exist
  
  checkFileExists(par::tpedfile);
  checkFileExists(par::tfamfile);


  ///////////////////////////////////////////////
  // .tfam file

  readFamFile(par::tfamfile);


  ///////////////////////////////////////////////
  // .tped file (map and genotype file)

  // Read through twice -- once as MAP (first four columns)
  // Then take genotypes
  
  // First rows columns of .tped file:

  // chromosome code
  // SNP identifier   
  // cM / M      
  // Base Position (-ve implies exclude)
  
  vector<bool> include;
  vector<Locus> ordered;
    
  FILE * MAP;
  MAP = fopen64(par::tpedfile.c_str(),"r");
  if ( MAP == NULL )
    error("Problem opening TPED file, errno = "+int2str(errno));

  int c=0;
  while( ! feof(MAP) )
    {
     
      // read first 4 entries, then ignore rest of line
      
      string chr;
      string name;
      string cm;
      string bp;

      long int inc;

      int f = 0;

      // Read chromosome code
      if ( readString(MAP,chr )) f++;
      
      // Empty line?
      if (chr=="") 
	continue;

      Locus * loc = new Locus;
      
      if ( readString(MAP,name )) f++;
      if ( readString(MAP,cm )) f++;
      if ( readString(MAP,bp )) f++;
            
      // Ignore rest of line
      while (fgetc(MAP) != '\n' && !feof(MAP)) {}

      // Store in locus information
      loc->name = name;
      loc->pos = atof(cm.c_str());
      loc->bp = (long int)atoi(bp.c_str());
      
      inc = loc->bp;

      // Check that cM/M specification looks correct, if 
      // we want to perform a plink-based analysis
      if (par::plink && (!par::cm_map) && (loc->pos > 50) )
	error("Looks like you need to specify --cm ??");
      
      // Convert cM to M map distances
      if (par::cm_map) loc->pos /= 100;
      
      // Chromosome coding
      loc->chr = getChromosomeCode(chr);

      // Use the frequency slot temporarily to 
      // store order information
      loc->freq = c++;
      
      // Are we including this locus?
      if (loc->name!="") 
	{
	  if (inc<0) 
	    {
	      include.push_back(false);
	    }
	  else 
	    {
	      include.push_back(true);
	      locus.push_back(loc);
	    }

	  ordered.push_back(*loc);
	}

      
    }
  
  // Done extracting initial SNP list
  fclose(MAP);
  
  printLOG(int2str(locus.size()) + " (of " + int2str(include.size()) + 
	   ") markers to be included from [ " + par::tpedfile + " ]\n");
  
  if ( locus.size() == 0 ) 
    shutdown();
  
  
  ///////////////////////////////////////////////
  // Build ordered table, so that genotypes can 
  // be inserted in correct order; then swap locus 
  // file over

  // Sorting a vector of pointers, so we need this special fix
  stable_sort(locus.begin(),locus.end(),less<Locus*>());
  
  // Sorting a normal vector
  stable_sort(ordered.begin(),ordered.end());


  c=0;
  for (int i=0; i<include.size(); i++)
    {
      // swap file order into locus position
      // and make all same chromosome
      ordered[i].bp = (int)ordered[i].freq;
      ordered[i].chr = 1;

      if (include[(int)ordered[i].bp]) 
	{
	  // keep track of genetic order, but only 
	  // for nonmissing loci
	  ordered[i].freq = c++;
	}      
      else ordered[i].freq = -1;
    }


  // resort to get lookup table
  stable_sort(ordered.begin(),ordered.end());  
  // i.e. for file position k, the locus position is ordered[k]->afreq


  // p2 p3 p1 p5 p4   : genetic position
  // 0  1  2  3  4    : file order
  // 1  0  1  0  1    : include

  // sort by cM
  // p1 p2 p3 p4 p5   : genetic 
  // 2  0  1  4  3    : file order
  // 1  1  0  1  0    : include
  // 0  1     2       : add genetic order: nonmissing...
  // 

  // sort by file order again
  // p2 p3 p1 p5 p4   : genetic
  // 0  1  2  3  4    : file
  // 1  0  1  0  1    : include
  // 1     0     2    : position to put in locus[l]
  
  

  ///////////////////////////////////////////////
  // Do we want to look at all the data?

  vector<int> include_pos(0);
  int nl_actual = locus.size();
  
  if ( (!par::plink) && (!par::run_chr==0) ) 
    {
      // Get range
      setMarkerRange();						
      
      // And set to 'exclude' all markers outside of this range 
      // (in physical distance terms)
      
      nl_actual = 0;
      
      for (int j=0; j<ordered.size(); j++)
	{
	  int fp = (int)ordered[j].freq;

	  if ( include[j] ) // not already excluded
	    {
	      // Is this SNP in range?
	      if ( fp < par::run_start || fp > par::run_end ) 
		{
		  include_pos.push_back(-1);
		  include[j] = false;
		}
	      else
		{
		  include_pos.push_back(fp);
		  nl_actual++;
		}
	    }
	  else // if already excluded
	    {
	      include_pos.push_back(-1);
	    }
	}
      
      //              0  1  2  3  4  5  6   7  8  9
      // We now have -1 -1 -1  3  4  5  6  -1 -1 -1
      // but we want -1 -1 -1  0  1  2  3  -1 -1 -1
      
      for (int j=0; j<ordered.size(); j++)
	{
	  if ( include_pos[j] > -1 ) 
	    include_pos[j] -= par::run_start ;
	}

    }
    else
    {
      // If we do want to look at all the data
      for (int j=0; j<ordered.size(); j++)
	include_pos.push_back((int)ordered[j].freq);
    }

  // Drop any markerspace from excluded range
  if ( (!par::plink) && (!par::run_chr==0) && nl_actual < locus.size() ) 
    {
      

      ////////////////////////////////////////
      // Remove selected loci from locus list, 
      // by copying rest to a new list
      
      vector<Locus*> l0(0);
      
      for(int l=0; l < locus.size(); l++)
	{    
	  
	  // If not in range	  
	  
	  if ( l < par::run_start || l > par::run_end ) 
	    {
	      // Free memory for original element
	      delete locus[l];
	    }
	  else
	    {
	      l0.push_back(locus[l]);
	    }
	}
      
      /////////////////
      // And copy back
      
      locus.clear();
      locus = l0;
      
    }


  ///////////////////////////////////////////////////
  // Add necessary locus space, if in SNP-major mode
  
  if (par::SNP_major)
    {
      for (int i=0; i<nl_actual; i++)
	{
	  CSNP * newlocus = new CSNP;
	  SNP.push_back(newlocus);
	}
    }
  else
    {
      for (int i=0; i<sample.size(); i++)
	{
	  sample[i]->one.resize(nl_actual);
	  sample[i]->two.resize(nl_actual);
	}
    }
  
  

  ///////////////////////////////////////////////
  // .tped, take 2 
  
  // Now re-read in the .tped file
  // Assume: the order of the FAM file specifies the
  // genotypes going across rows, so we should already 
  // have that information

  FILE * PED;
  PED = fopen64(par::tpedfile.c_str(),"r");
  
  int i=0; // SNP count

  while( ! feof(PED) )
    {
      
      string dummy;
      int f=0;

      if (readString(PED,dummy )) f++;
      
      // End of file?
      if ( dummy=="" )
	{
	  continue;
	}
      
      // Is this line a comment, or are we skipping it?
      if ( dummy.substr(0,1)=="#" )
	{
	  // Ignore rest of line
	  while (fgetc(PED) != '\n' && !feof(PED)) {}	  
	  continue;
	}


      if ( !include[i] )
	{
	  // Ignore rest of line and advance to next SNP
	  while (fgetc(PED) != '\n' && !feof(PED)) {}
	  i++;
	  continue;
	}



      // Skip next 3 fields (SNP, cM, bp)

      if (readString(PED,dummy)) f++;
      if (readString(PED,dummy)) f++;
      if (readString(PED,dummy)) f++;
            
      
      /////////////////////
      // Read genotypes now
      
      int gn=0;
      int c=0; // individual count
      bool linedone = false;
      bool fatal = false;
      string fmsg;
      while ( ! linedone )
	{
	  
	  Individual * person = sample[c];
	  
	  string one="";
	  string two="";
	  
	  while (1)
	    {
	      char ch = fgetc(PED);
	   
	      // Delimiter?
	      if (ch==' ' || ch=='\t' || ch=='\n' || ch=='\r' || feof(PED) )
		{
		  
		  if (ch=='\n' || ch=='\r' || feof(PED))
		    linedone = true;
		  
		  // have we already seen something?
		  if (one.length()>0) 
		    {
		      gn++;		  
		      break;
		    }
		  
		  if (ch=='\n' || ch=='\r' || feof(PED))
		    break;
		  
		}
	      else
		{
		  one += ch;
		}
	    }
	  
	  // Second allele
	  if (!linedone)
	    while (1)
	      {
		char ch = fgetc(PED);
		
		// Delimiter?
		if (ch==' ' || ch=='\t' || ch=='\n' || ch=='\r' || feof(PED) )
		  {
		    
		    if (ch=='\n' || ch=='\r' || feof(PED))
		      linedone = true;
		    
		    // have we already seen something?
		    if (two.length()>0) 
		      {
			gn++;      
			break;
		      }
		    
		    if (ch=='\n' || ch=='\r' || feof(PED))
		      break;
		    
		  }
		
		else
		  {
		    two += ch;
		  }
	      }
	  
	  
	  if (linedone && one.length()==0 && two.length()==0 )
	    break;
	  

	  /////////////////////////////////////
	  // Only consider loci to be included
	  
	  if (include[i]) 
	    {
	      
	      //////////////////////////////
	      // Look up genomic order, 
	      // insert in slot k in locus[]

	      int k = include_pos[i];

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

	      
	      /////////////////////////////
	      // Add specific genotypes
	      
	      if (par::SNP_major)
		{
		  
		  // 00 hom
		  if (one==loc->allele1 && two==loc->allele1)
		    {
		      SNP[k]->one.push_back(false);
		      SNP[k]->two.push_back(false);
		    }
		  
		  
		  // 01 het
		  else if (one!=par::missing_genotype && 
			   two!=par::missing_genotype && 
			   one!=two)	    
		    {
		      SNP[k]->one.push_back(false);
		      SNP[k]->two.push_back(true);
		    }
		  
		  // 11 hom
		  else if (one==loc->allele2 && two==loc->allele2)
		    {
		      SNP[k]->one.push_back(true);
		      SNP[k]->two.push_back(true);
		    }
		  
		  // 10 missing
		  else if (one==par::missing_genotype || two==par::missing_genotype)
		    {
		      SNP[k]->one.push_back(true);
		      SNP[k]->two.push_back(false);
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


	  // Advance to next individual
	  c++;
	  
	  if ( c > sample.size())
	    {
	      fmsg += "\nProblem with line "+int2str(i+1)+" in [ "+par::tpedfile+" ]\n";
	      fmsg += "Expecting 4 + 2 * " + int2str(sample.size()) + " = " + 
		int2str(4+2*sample.size())+ " columns, but found more\n";
	      error(fmsg);
	    }
	  
	} // line done? Next SNP
      
      
      // check size of line length somewhere
      if ( gn != 2 * sample.size() )
	{	  
	  fmsg += "\nProblem with line "+int2str(i+1)+" in [ "+par::tpedfile+" ]\n";
	  fmsg += "Expecting 4 + 2 * " + int2str(sample.size()) + " = " + 
	    int2str(4+2*sample.size())+ " columns, but found " + 
	    int2str(f+gn) + "\n";
	  fatal=true;
	}
      

      if (fatal) 
	error(fmsg);

      // Increase SNP counter
      i++;
      

    } // Next SNP
  
  // Close TPED file
  fclose(PED);


}  



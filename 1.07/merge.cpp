

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




void Plink::mergeList()
{
  
  // Perform multiple merge operations (first fileset is specified on
  // the command line, the rest in a merge-list). Simply perform
  // multiple calls to mergeData()

  // If a line contains 2 items, it suggests a PED/MAP merge
  // If a line contains 3 items, it suggests a BED/BIM/FAM merge

  printLOG( "Using merge mode " + int2str( par::merge_mode ) + " : ");
  if (par::merge_mode==1) printLOG("consensus call (default)\n");
  else if (par::merge_mode==2) printLOG("overwrite if missing in original\n");
  else if (par::merge_mode==3) printLOG("overwrite unless missing in new\n");
  else if (par::merge_mode==4) printLOG("overwrite none\n"); 
  else if (par::merge_mode==5) printLOG("overwrite all\n");

  
  // Check that a merge-list exists
  checkFileExists(par::merge_list_filename);
  ifstream MLIST(par::merge_list_filename.c_str());
  MLIST.clear();


  // Iterate over all files to be merged
  int c=0;
  while(!MLIST.eof())
    {

      char cline[5000] = "";
      MLIST.getline(cline,5000,'\n');

      // convert to string
      string sline = cline;
      if (sline=="") continue;
      
      string buf; 
      stringstream ss(sline); 
      vector<string> tokens; 
      while (ss >> buf)
	tokens.push_back(buf);

	// Is this a PED/MAP, or a BED/BIM/FAM merge?      
      if (tokens.size() == 2) 
	{
	  par::merge_pedfile = tokens[0];
	  par::merge_mapfile = tokens[1];
	  if (par::merge_pedfile=="" || par::merge_mapfile=="") continue;
	  par::merge_binary = false;
	} 
      else if (tokens.size() == 3)
	{
	  par::merge_bedfile = tokens[0];
	  par::merge_bimfile = tokens[1];
	  par::merge_famfile = tokens[2];
	  if (par::merge_bedfile=="" || par::merge_bimfile=="" || par::merge_famfile=="") continue;      
	  par::merge_binary = true;
	}
      else error("Problem with merge-list file line:\n should be either 2 (PED/MAP) or 3(BED/BIM/FAM) fields:\n"+sline);
      
      c++;
      if (!par::silent)
	cout << c << " files merged            \r";
      
      // Perform the actual merge
      if (!par::merge_binary) mergeData();
      else mergeBinaryData();
      
      // Reset number of individuals
      n = sample.size();
      
      // Set number of pairs
      np = (int)((double)(n*(n-1))/(double)2);
      
      // Total number of all (test+background) loci
      nl_all = locus.size();
      
    }
  
  if (!par::silent)
    cout << "\n";
  MLIST.close();
  
  // Report some stats now we've finished
  printLOG("Merging "+int2str(c)+" samples, final sample contains "+int2str( n ) + " individuals");
  printLOG(" and " +int2str( nl_all ) + " markers\n");
  
}



void Plink::mergeData()
{
  
  if (!par::SNP_major)
    error("--merge can only apply in SNP-major mode currently\n");

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
    }
  
  // We've already loaded in the first file
  // Do not overwrite any existing phenotype information
  
  checkFileExists(par::merge_pedfile);
  checkFileExists(par::merge_mapfile);
  
  // Make hash of original SNP names
  
  map<string,int> mlocus;
  for (int l=0;l<nl_all;l++)
    mlocus.insert(make_pair(locus[l]->name,l));
  map<string,int>::iterator ilocus;
  

  // Reset counts
  diff_overlap = 0;
  diff_nonmissing_overlap = 0;
  diff_concordant_overlap = 0;



  // A temporary hash for the names of any markers that
  // do not match in terms of strand
  map<string,int> misstrand;


  ///////////////////////////////////////
  // .map 

  vector<bool> include(0);
  map<string,bool> exists;
  vector<Locus> ordered(0);
  vector<Locus*> locus2(0);
  
  ifstream MAP(par::merge_mapfile.c_str());
  MAP.clear();

  int exist_cnt=0;
  
  int c=0;
  while(!MAP.eof())
    {
      
      Locus * loc = new Locus;
      
      string chr;
      long int inc;
      
      char cline[256] = "";
      MAP.getline(cline,256,'\n');

      // convert to string
      string sline = cline;
      if (sline=="") continue;

      string buf; 
      stringstream ss(sline); 
      vector<string> tokens; 
      while (ss >> buf)
	tokens.push_back(buf);
      
      if (tokens.size() == 0)
	continue;
      else if ( par::map3 && tokens.size() != 3 ) 
	error("Problem with MAP file line:\n"+sline);
      else if ( (!par::map3) && tokens.size() != 4 )
	error("Problem with MAP file line:\n"+sline);
      
      chr = tokens[0];
      loc->name = tokens[1];
      
      if ( par::map3 )
	{
	  loc->pos = 0;
	  loc->bp = (long int)atoi(tokens[2].c_str());
	}
      else
	{
	  loc->pos = atof(tokens[2].c_str());
	  loc->bp = (long int)atoi(tokens[3].c_str());
	}
      
      inc = loc->bp;
      
      // Check that cM/M specification looks correct, if 
      // we want to perform a plink-based analysis
      if (par::plink && !par::cm_map && loc->pos > 50)
	error("Looks like you need to specify --cm ??");
      
      // Convert cM to M map distances
      if (par::cm_map) loc->pos /= 100;
      
      // Convert chromosome code, taking species into account
      loc->chr = getChromosomeCode( chr );

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
	      ilocus = mlocus.find(loc->name);
	      
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
		  include.push_back(true);
		  locus2.push_back(loc2);
		  
		  // Keep the new file order (would have been in freq)
		  int t = (int)loc->freq;
		  
		  // Clean up if we do not need another locus
		  delete loc;
		  
		  // Replace
		  loc = loc2;
		  loc->freq = t;
		}
	      else 
		{
		  // Locus does not exist -- add to locus list
		  exists[loc->name] = false;
		  include.push_back(true);
		  locus2.push_back(loc);
		}
	      
	    }
	  
	  ordered.push_back(*loc);
	}
      
      
    }

  MAP.clear();
  MAP.close();
  
  if (!par::merge_list)
    {
      printLOG("\n" +int2str(locus2.size()) + 
	       " (of " + int2str(include.size()) +
	       ") markers to be merged from [ " +par::merge_mapfile + " ]\n");
      printLOG("Of these, "+int2str( locus2.size()-exist_cnt ) + " are new, " +
	       int2str( exist_cnt ) +  " already exist\n");
    }

  ///////////////////////////////////////////////
  // Build ordered table, so that genotypes can be inserted 
  // in correct order; then swap locus file over
  
  // Sorting a vector of pointers, so we need this special fix
  stable_sort(locus2.begin(),locus2.end(),less<Locus*>());
  
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
      else {
	ordered[i].freq = -1;
      }

    }
  
  // resort to get lookup table
  stable_sort(ordered.begin(),ordered.end());  

  // i.e. for file position k, the locus position is ordered[k]->freq

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
  // .ped
  
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
  // Read in PED file to new merge file
  
  bool fatal_error = false;

  FILE * PED;
  PED = fopen(par::merge_pedfile.c_str(),"r");

  c=0;
  while(! feof(PED) )
    {
      
      Individual * person = new Individual;

      // Get first field
      int f=0;
      if (readString(PED,person->fid )) f++;
      
      // End of file?
      if ( person->fid=="" )
	continue;

      // Check for reserved family ID code
      if ( person->fid=="FID" )
	error("FID is a reserved ID... please select a different family ID");

      // Is this line a comment?
      
      if (person->fid.substr(0,1)=="#")
	{
	  // Ignore rest of line and advance to next line
	  while (fgetc(PED) != '\n' && !feof(PED)) {}
	  continue;
	}
      
      // First 6 obligatory fields

      if ( readString(PED,person->iid )) f++;
      if ( readString(PED, person->pat )) f++;
      if ( readString(PED, person->mat )) f++;
      if ( readString(PED, person->sexcode)) f++;

      string phenotype;
      if (readString(PED,phenotype)) f++;

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


      // Skip last empty line that gets read
      if (person->fid=="") break;


      if (person->sexcode=="1")
	person->sex = true;  // male
      else if (person->sexcode=="2")
	person->sex = false; // female (default)

      
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
	      if ( ! from_string<double>( person->phenotype, phenotype, std::dec ))
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
	  
// 	  if (!par::merge_list)
// 	    if (person->missing) 
// 	      {
// 		stringstream s2;
// 		s2 << "Individual " << person->fid << " " 
// 		   << person->iid << " has missing phenotype: "
// 		   << person->phenotype << " / " 
// 		   << par::missing_phenotype << "\n";
// 		printLOG(s2.str());
// 	      }
	  
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

      

      
      /////////////////////////
      // For each locus

      int gn=0;
      int i=0;
      bool linedone = false;
      bool fatal = false;
      string fmsg;
      while ( ! linedone )
	{

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

	      int k0 = (int)ordered[i].freq;
	      ilocus = mlocus.find(locus2[k0]->name);
	      int k = ilocus->second;

	      bool e = reconcileMerge(indn,
				      k,
				      one,two,
				      already_in, 
				      exists.find(locus2[k0]->name)->second,
				      MERD,
				      misstrand );
	      
	      if (e) fatal_error = true;
	      
	    }
	  
	  // Advance to next locus
	  i++;
	  if ( i > include.size())
	    {
	      fmsg += "\nProblem with line "+int2str(c+1)+" in [ "
		+par::merge_pedfile+" ]\n";
	      fmsg += "Expecting 6 + 2 * " + int2str(include.size()) + " = " + 
		int2str(6+2*include.size())+ " columns, but found more\n";
	      error(fmsg);
	    }	  
	} // Next locus
      
      // check size of line length somewhere
      if ( gn != 2 * include.size() )
	{	  
	  fmsg += "\nProblem with line "+int2str(c+1)+" in [ "
	    +par::merge_pedfile+" ]\n";
	  fmsg += "Expecting 6 + 2 * " + int2str(include.size()) + " = " + 
	    int2str(6+2*include.size())+ " columns, but found " + 
	    int2str(f+gn) + "\n";
	  fatal=true;
	}
      
      if (fatal) 
	error(fmsg);
      
      // Increase person counter
      c++;
      
      // Add individual to list, if need be
      if (!already_in)
	sample.push_back(person);
      
    }

  
  ///////////////////////////////////////////////////////
  // Did we encounter any fatal errors from flipped SNPs? 
  
  if (fatal_error)
    {
      
      ofstream MSNP;
      string f = par::output_file_name+".missnp";
      MSNP.open(f.c_str(), ios::out);
      map<string,int>::iterator ilocus;
      
      for ( ilocus = misstrand.begin() ; ilocus != misstrand.end() ; ilocus++)
	{
	  MSNP << locus[ilocus->second]->chr << "\t"
	       << locus[ilocus->second]->name << "\n";
	}
      MSNP.close();
      cerr << "\nFound " << misstrand.size() << " SNPs that do not match in terms of allele codes\n";
      cerr << "Might include strand flips, although flipped A/T and C/G SNPs will be undetected)\n";
      cerr << "Writing problem SNPs to [ " << f << " ]\n";
      error("Stopping due to mis-matching SNPs -- check +/- strand?");
    }
  
  
  // If a binary trait, now make 0 missing also
  // i.e. if we never saw other than missing, 0, 1 or 2
  
  if (par::bt)
    for (int i=0; i<sample.size(); i++)
      if ( sample[i]->phenotype == 0 )
	sample[i]->missing = true;


  // Close the PED file
  fclose(PED);

  
  if (!par::merge_list)
    {
      printLOG(int2str( c ) + " individuals merged from [ " 
	       + par::merge_pedfile + " ] \n");
      printLOG("Of these, " + int2str(new_person) + " were new, " 
	       + int2str( old_person ) +
	       " were already in current data\n\n");
    }
  
  if (par::merge_mode >=6) 
    {

      printLOG("Results from diff ( merge mode " 
	       + int2str(par::merge_mode) + " ) written to [ " +
	       par::output_file_name + ".diff ]\n");

      MERD.close();

      printLOG("Of " + int2str(diff_overlap) 
	       + " overlapping SNPs, " 
	       + int2str( diff_nonmissing_overlap ) 
	       + " were both genotyped\nand "
	       + int2str( diff_concordant_overlap ) 
	       + " were concordant\n");
      printLOG("Concordance rate is " 
	       + dbl2str( (double)diff_concordant_overlap / (double)diff_nonmissing_overlap ) + "\n");
      
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





bool Plink::reconcileMerge(int indn, int k, string one, string two, 
			   bool already_in, 
			   bool snp_exists, 
			   ofstream & MERD,
			   map<string,int> & misstrand)
{

//   cout << "rec locus " << locus[k]->name << "\n";
//   cout << "in = " << one << " " << two <<"\n";
//   cout << "existing = " << locus[k]->allele1 
//        << " " << locus[k]->allele2 << "\n";
  
  // Note -- this routine does not fully implemented individual-major
  // mergeing, which is why we currently only allow SNP-major in the 
  // main mergeData() routine above
  
  bool fatal_error = false;
  
  Locus * loc = locus[k];
  Individual * person = NULL;
  if (already_in) person = sample[indn];

  
  /////////////////////////////////////////
  // Add allele names to list, if needed
  // 
  // For a merged-in binary file, we already will have performed
  // this step when reading in the BIM file

  if ( ! par::merge_binary ) 
    {
      
      // If allele is not missing...
      if (one!=par::missing_genotype && two!=par::missing_genotype)
	{
	  // ...and not already listed
	  if (one!=loc->allele1 && one!=loc->allele2)
	    {
	      // ...then add to first empty slot
	      if(loc->allele1=="" || loc->allele1==par::missing_genotype)
		loc->allele1=one;
	      else if(loc->allele2=="" || loc->allele2==par::missing_genotype)
		loc->allele2=one;
	      else {
		// .. or show an error if no empty slots
		misstrand.insert(make_pair(loc->name,k));
		fatal_error = true;
	      }
	    }
	}
      
      
      //////////////////////////////////////////
      // Repeat for second allele, if different
      
      if (two!=one)
	{
	  // If allele is not missing...
	  if (two!=par::missing_genotype)
	    // ...and not already listed
	    if (two!=loc->allele1 && two!=loc->allele2)
	      {
		// ...then add to first empty slot
		if(loc->allele1=="" || loc->allele1==par::missing_genotype) 
		  loc->allele1=two;
		else if(loc->allele2=="" || loc->allele2==par::missing_genotype ) 
		  loc->allele2=two;
		else 
	      {
		misstrand.insert(make_pair(loc->name,k));
		fatal_error = true; 
	      }
	      }
	}
      
    }
  
  /////////////////////////
  // Add specific genotypes
  
  bool write;
  
  if ( !already_in ) write = true;        // Write if new person
  else if ( ! snp_exists ) write = true;  // Write if new SNP
  else 
    {

      // The genotype exists in both the original and the new files
      // Depending on the merge-mode, we need to determine whether
      // to overwrite, or report, these genotypes
      
      bool s1;
      bool s2;
      
      if (par::SNP_major)
	{
	  s1 = SNP[k]->one[indn];
	  s2 = SNP[k]->two[indn];
	}
      else
	{
	  s1 = person->one[k];
	  s2 = person->two[k];
	}
      
      // MODE 1: Consensus call
      if ( par::merge_mode == 1) 
	{
	  // If new genotype missing, never write
	  if ( one==par::missing_genotype || two==par::missing_genotype) write = false;
	  
	  // If existing is missing, always write
	  else if ( (s1) && (!s2) ) write = true;
	  
	  // Else if both called, check they match
	  else
	    {
	      
	      bool mismatch = false;
	      if ( one==loc->allele1 && two==loc->allele1 ) // New == 11
		{
		  if ( ! ( (!s1) && (!s2) ) ) mismatch = true;
		}
	      else if ( one==loc->allele1 && two==loc->allele2 ) // New == 12
		{
		  if ( ! ( (!s1) && s2 ) ) mismatch = true;
		}
	      else if ( one==loc->allele2 && two==loc->allele1 ) // New == 12
		{
		  if ( ! ( (!s1) && s2 ) ) mismatch = true;
		}
	      else if ( one==loc->allele2 && two==loc->allele2 ) // New == 22
		{
		  if ( ! ( (s1) && (s2) ) ) mismatch = true;
		}
	      
	      if (mismatch) 
		{
		  one = par::missing_genotype;
		  two = par::missing_genotype;
		  write = true;
		}
	    }
	}  
      
      // MODE 2: Overwrite if original missing
      else if ( par::merge_mode == 2) 
	{
	  if ( s1 && (!s2) ) write = true;
	  else write = false;
	}
      
      // MODE 3: Overwrite unless missing in new 
      else if ( par::merge_mode == 3)
	{
	  if ( one==par::missing_genotype || two==par::missing_genotype) write = false;
	  else write = true;
	}
      
      // MODE 4: Never overwrite
      if ( par::merge_mode == 4) write = false;
      
      // MODE 5: Overwrite all
      else if ( par::merge_mode == 5)
	write = true;
      
      // MODE 6,7 : Report diffs (if non-missing)
      else if ( par::merge_mode == 6 || par::merge_mode == 7 )
	{
	  
	  bool mismatch = false;
	  
	  bool new_geno_missing = one==par::missing_genotype || two==par::missing_genotype ;
	  
	  if ( new_geno_missing  ) // New == ??
	    {
	      if ( ! ( s1 && (!s2) ) ) mismatch = true;
	    }
	  else if ( one==loc->allele1 && two==loc->allele1 ) // New == 11
	    {
	      if ( ! ( (!s1) && (!s2) ) ) mismatch = true;
	    }
	  else if ( one==loc->allele1 && two==loc->allele2 ) // New == 12
	    {
	      if ( ! ( (!s1) && s2 ) ) mismatch = true;
			}
	  else if ( one==loc->allele2 && two==loc->allele1 ) // New == 12
	    {
	      if ( ! ( (!s1) && s2 ) ) mismatch = true;
	    }
	  else if ( one==loc->allele2 && two==loc->allele2 ) // New == 22
	    {
	      if ( ! ( (s1) && (s2) ) ) mismatch = true;
	    }
	  
	  if (par::merge_mode==7) 
	    {
	      if ( new_geno_missing || ( s1 && (!s2) ) )
		mismatch=false;	      
	    }

	  // Summary stats
	  ++diff_overlap;
	  if ( ! ( new_geno_missing || ( s1 && (!s2) ) ) )
	    {
	      ++diff_nonmissing_overlap;
	      if ( ! mismatch ) 
		++diff_concordant_overlap;
	    }
	  
	  if (mismatch)
	    {
	      MERD << setw(20) << loc->name << " "
		   << setw(20) << person->fid << " " 
		   << setw(20) << person->iid << " "
		   << setw(8) << (string)(one+"/"+two) << " ";
	      
	      if ((!s1) && (!s2))
		MERD << setw(8) << (string)(loc->allele1+"/"+loc->allele1) << " ";
	      
	      if ((!s1) && s2)
		MERD << setw(8) << (string)(loc->allele1+"/"+loc->allele2) << " ";
	      
	      if ( s1 && s2 )
		MERD << setw(8) << (string)(loc->allele2+"/"+loc->allele2) << " ";
	      
	      if ( s1 && (!s2))
		MERD << setw(8) << (string)(par::missing_genotype+"/"+par::missing_genotype) << " ";
	      
	      MERD << "\n";
	      
	    }
	  

	}
      
    }
  
	     

  if (write)
    {

      
      if (par::SNP_major)
	{

	  // 00 hom
	  if (one==loc->allele1 && 
	      two==loc->allele1 &&
	      one!=par::missing_genotype )
	    {
	      SNP[k]->one[indn] = false;
	      SNP[k]->two[indn] = false;
	    }
		      
	  
	  // 01 het
	  else if (one!=par::missing_genotype && 
		   two!=par::missing_genotype && 
		   one!=two)	    
	    {
	      SNP[k]->one[indn] = false;
	      SNP[k]->two[indn] = true;
	    }
	  
	  // 11 hom
	  else if (one==loc->allele2 && 
		   two==loc->allele2 &&
		   one!=par::missing_genotype)
	    {
	      SNP[k]->one[indn] = true;
	      SNP[k]->two[indn] = true;
	    }
	  
	  // 10 missing
	  else 
	    {
	      SNP[k]->one[indn] = true;
	      SNP[k]->two[indn] = false;
	    }
	}
      else
	{
	  /////////////////////////
	  // Individual-major mode
	  
	  // 00 hom
	  if (one==loc->allele1 && 
	      two==loc->allele1 && 
	      one!=par::missing_genotype)
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
	  else if (one==loc->allele2 && 
		   two==loc->allele2 &&
		   one!=par::missing_genotype)		   
	    {
	      person->one[k]=true;
	      person->two[k]=true;		  
	    }
	  
	  // 10 missing
	  else
	    {
	      person->one[k]=true;
	      person->two[k]=false;		  
	    }
	}
    }	      
  

  return fatal_error;
}



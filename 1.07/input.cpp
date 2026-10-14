

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
#include <cerrno>

#include "plink.h"
#include "options.h"
#include "helper.h"
#include "nlist.h"
#include "gvar.h"

extern ofstream LOG;

void Plink::readData()
{
  
  //////////////////////
  // Check files exist
  
  if ( ! par::ped_from_stdin)
    checkFileExists(par::pedfile);

  checkFileExists(par::mapfile);

  ///////////////////////////////////////////////
  // .map file

  vector<bool> include;
  vector<int> include_pos(0);
  int nl_actual=0;
  
  readMapFile(par::mapfile,
	      include,
	      include_pos,
	      nl_actual);


  
  ///////////////////////////////////////////////
  // .ped
  

  FILE * PED;

  if ( ! par::ped_from_stdin )
    {
      PED = fopen64(par::pedfile.c_str(),"r");
      if ( PED == NULL )
	error("Problem opening PED file, errno = "+int2str(errno));
    }

  vector<Individual*> ambiguous;

  int nmale = 0;
  int nfemale = 0;
  int nambig = 0;

  int c=0; // number of individuals
  string s2;
  

  while( 1 )
    {

      // End of input stream? 

      if ( par::ped_from_stdin )
	{
	  if ( cin.eof() )
	    break;
	}
      else
	{
	  if ( feof(PED) )
	    break;
	}


      // Otherwise read in the next person

      Individual * person = new Individual;
      

      // Get first field
      int f=0;
      
      if ( par::ped_from_stdin)
	cin >> person->fid;
      else
	if (readString(PED,person->fid )) f++;
      
      // End of file?
      if ( person->fid=="" )
	{
	  delete person;
	  continue;
	}

      if ( person->fid=="FID" )
	error("FID is a reserved ID... please select a different family ID");
      

      // Is this line a comment?      
      if ( ! par::ped_from_stdin)
	{
	  if (person->fid.substr(0,1)=="#")
	    {
	      // Ignore rest of line and advance to next line
	      while (fgetc(PED) != '\n' && !feof(PED)) {}
	      delete person;
	      continue;
	    }
	}


	// First 6 or 7 obligatory fields
      
      if ( par::ped_skip_fid )
	person->iid = person->fid;
      else
	{
	  if ( par::ped_from_stdin )
	    cin >> person->iid;
	  else
	    if ( readString(PED,person->iid )) f++;
	}

      if ( par::ped_skip_parents )
	{
	  person->mat = person->pat = "0";
	}
      else
	{
	  if ( par::ped_from_stdin )
	    {
	      cin >> person->pat >> person->mat;
	    }
	  else
	    {
	      if ( readString(PED, person->pat )) f++;
	      if ( readString(PED, person->mat )) f++;
	    }
	}
      
      if ( par::ped_skip_sex )
	person->sexcode = "0";
      else
	{
	  if ( par::ped_from_stdin )
	    cin >> person->sexcode;
	  else
	    if ( readString(PED, person->sexcode)) f++;
	}
      
      string phenotype;
      if ( par::ped_skip_pheno )
	phenotype = par::missing_phenotype;
      else
	{
	  if ( par::ped_from_stdin )
	    cin >> phenotype;
	  else
	    if (readString(PED,phenotype)) f++;
	}



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
      
      // Optional liability class
      if (par::liability)
 	{
 	  string dummy;
	  if ( par::ped_from_stdin )
	    cin >> dummy;
	  else
	    if (readString(PED,dummy)) f++;	    
 	}

      
      // Skip last empty line that gets read
      if (person->fid=="") break;

      // Check sex
      if (person->sexcode=="1")
	{
	  person->sex = true; // male
	  nmale++;
	}
      else if (person->sexcode=="2")
	{
	  person->sex = false;  // female
	  nfemale++;
	}
      else 
	{
	  ambiguous.push_back(person);
	  nambig++;
	  if (!par::ignore_missing_sex)
	    person->missing = true;
	}
      

      //////////////////
      // A non-founder?
      
      person->founder = 
	(person->pat == "0" && person->mat == "0") ? true : false;
      

      //////////////////////////////
      // Test for quantitative trait
      
      if (phenotype == par::missing_phenotype)
	person->missing = true;
      else
	{
	  
	  // Store in person->phenotype as number, checking for 
	  // conversion failure
	  if (  ! from_string<double>( person->phenotype, phenotype, std::dec ) )
	    person->missing = true;
	  else
	    {
	      if (phenotype != "0" && 
		  phenotype != "1" && 
		  phenotype != "2" ) 
		{
		  par::qt = true;
		  par::bt = false; 
		}
	    }
	}
      
            
      /////////////////////////////
      // Add necessary locus space

      if (!par::SNP_major)
	{
	  person->one.resize(nl_actual);
	  person->two.resize(nl_actual);
	}
      
      /////////////////////
      // Read genotypes now
      
      int gn=0;
      int i=0;
      bool linedone = false;
      bool fatal = false;
      string fmsg;
      while ( ! linedone )
	{
	  
	  string one="";
	  string two="";
	
	  if ( par::ped_from_stdin )
	    cin >> one >> two;
	  else
	    {
  
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
	      
	      // Is this a compound genotype?
	      if ( par::compound_genotype_code )
		{
		  // In this case, each allele must be exactly 1-character long
		  if ( one.length() != 2 ) 
		    error("Problem with compound genotype [ " + one + " ] should be two characters long\n");
		  
		  two = one[1];
		  one = one[0];

		  // Add second allele
		  ++gn;
		}
	      else
		{
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
		}
	      
	      if (linedone && one.length()==0 && two.length()==0 )
		break;
	      
	    }
	  
	  
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
			    + person->fid + " " + person->iid 
			    + " has genotype [ " + one +" "+two+" ]\n"
			    + "       but we've already seen [ " 
			    + loc->allele1 + " ] and [ " + loc->allele2 + " ]\n";
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
			    + person->fid + " " + person->iid 
			    + " has genotype [ " + one +" "+two+" ]\n"
			    + "       but we've already seen [ " 
			    + loc->allele1 + " ] and [ " + loc->allele2 + " ]\n";
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


	  // Advance to next locus
	  i++;

	  if ( par::ped_from_stdin )
	    {
	      if ( i == include.size() )
		linedone = true;
	    }
	  else
	    {
	      if ( i > include.size())
		{
		  int ef = 1;
		  if ( ! par::ped_skip_fid ) ef++;
		  if ( ! par::ped_skip_parents ) ef+=2;
		  if ( ! par::ped_skip_sex ) ef++;
		  if ( ! par::ped_skip_pheno ) ef++;
		  if (   par::liability ) ef++;
		  
		  fmsg += "\nProblem with line "
		    +int2str(c+1)+" in [ "+par::pedfile+" ]\n";
		  fmsg += "Expecting "
		    +int2str(ef) + " + 2 * " + int2str(include.size()) + " = " + 
		    int2str(ef+2*include.size())+ " columns, but found more\n";
		  error(fmsg);
		}
	    }

	} // line done?
    
      
      // check size of line length somewhere
      if ( ! par::ped_from_stdin )
	{
	  if ( gn != 2 * include.size() )
	    {	  
	      int ef = 1;
	      if ( ! par::ped_skip_fid ) ef++;
	      if ( ! par::ped_skip_parents ) ef+=2;
	      if ( ! par::ped_skip_sex ) ef++;
	      if ( ! par::ped_skip_pheno ) ef++;
	      if (   par::liability ) ef++;
	      
	      fmsg += "\nA problem with line "
		+int2str(c+1)+" in [ "+par::pedfile+" ]\n";
	      fmsg += "Expecting "+int2str(ef)
		+" + 2 * " + int2str(include.size()) + " = " + 
		int2str(ef+2*include.size())+ " columns, but found " + 
		int2str(f+gn) + "\n";
	      fatal=true;
	    }
	}

      if (fatal) 
	error(fmsg);

      // Increase person counter
      c++;
      
      // Add individual to list
      sample.push_back(person);
      
    }
  
  // If a binary trait, now make 0 missing also
  // i.e. if we never saw other than missing, 0, 1 or 2
  
  if (par::bt)
    for (int i=0; i<sample.size(); i++)
      if ( sample[i]->phenotype == 0 )
	sample[i]->missing = true;

  // Display list of ambiguously-sexed individuals?
  if (ambiguous.size()>0)
    {
      printLOG("Warning, found "
	       +int2str(ambiguous.size())
	       +" individuals with ambiguous sex codes\n");
      if (!par::ignore_missing_sex)
	printLOG("These individuals will be set to missing ( or use --allow-no-sex )\n");      
      string f = par::output_file_name + ".nosex";
      printLOG("Writing list of these individuals to [ "+f+" ]\n");
      ofstream AMB;
      AMB.open(f.c_str(), ifstream::out);
      for (int i=0; i<ambiguous.size(); i++)
	AMB << ambiguous[i]->fid << "\t" << ambiguous[i]->iid << "\n";
      AMB.close();      
      ambiguous.clear();
    }

     

  // Close PED file
  if ( ! par::ped_from_stdin )
    fclose(PED);

    
  
  printLOG(int2str(c)+" individuals read from [ "+par::pedfile+" ] \n");
  int nm=0;
  for (int i=0;i<sample.size();i++)
    if(!sample[i]->missing) nm++;
  printLOG(int2str(nm)+" individuals with nonmissing phenotypes\n");

  if (par::bt) 
    {

      if (par::coding01) 
	printLOG("Assuming a disease phenotype (0=unaff, 1=aff, other=miss)\n");
      else
	{
	  printLOG("Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)\n");
	  if (par::missing_phenotype!="0")
	    printLOG("Missing phenotype value is also " + par::missing_phenotype + "\n");
	}

      int ncase = 0;
      int ncontrol = 0;
      int nmissing = 0;

      for (int i=0; i<sample.size(); i++)
	if ( sample[i]->missing ) 
	  nmissing++;
	else if ( sample[i]->phenotype == 1 )
	  ncontrol++;
	else if ( sample[i]->phenotype == 2 )
	  ncase++;

      printLOG(int2str(ncase)+" cases, "
	       +int2str(ncontrol)+" controls and "
	       +int2str(nmissing)+" missing\n");
    }
  else 
    {
      printLOG("Assuming a quantitative trait\n");
      printLOG("Missing phenotype value is " 
	       + par::missing_phenotype + "\n");
    }

  // Display sex counts
  printLOG(int2str(nmale)+" males, "+int2str(nfemale)
	   +" females, and "+int2str(nambig)+" of unspecified sex\n");
  
  
}



void Plink::readSet()
{

  bool firsttime = true;

  set<string> subset;

  if ( par::use_subset )
    {
      printLOG("Reading a list of subsets from [ " + par::subsetfile + " ]\n");

      checkFileExists( par::subsetfile );
      ifstream IN(par::subsetfile.c_str(), ios::in);
  
      while ( ! IN.eof() )
	{
	  string gname;
	  IN >> gname;
	  if ( gname=="" )
	    continue;
	  subset.insert(gname);
	}
      printLOG("Read " + int2str( subset.size() ) + " sets to extract\n");
    }
  
      

  //////////////////////
  // Clear current sets

  // (i.e. after removing SNPs, the lookup numbers will be wrong so we
  // need to reload)
  
  if (snpset.size()>0)
    {
      firsttime = false;
      for (int i=0; i<snpset.size();i++)
	snpset[i].clear();
      snpset.clear();
      setname.clear();
  }
  
  //////////////////
  // Load SET file
  
  checkFileExists(par::setfile);
  ifstream SET;
  SET.open(par::setfile.c_str());
  SET.clear();
  
  // Temporary vector of SNP numbers
  vector<int> s;
  
  // First set name
  string name;
  SET >> name;
  
  // Make map of locus name with 'l' number
  map<string,int> mlocus;
  for (int l=0;l<nl_all;l++)
    mlocus.insert(make_pair(locus[l]->name,l));
  
  map<string,int>::iterator ilocus;
  
  
  while(!SET.eof())
    {
      string t;
      SET >> t;
      
      if (t=="END" || t=="end" || SET.fail() ) // End of SET
	{
	  
	  if ( SET.fail() ) 
	    printLOG("Warning: the set-file did not end with the END keyword\n"); 
	  
	  if ( ( ! par::use_subset ) ||
	       subset.find(name) != subset.end() )
	    {	
	      
	      // Save set
	      snpset.push_back(s);	  

	      // Save set name
	      setname.push_back(name);
	    }
	  
	  
	  // Get next set name
	  SET >> name;

	  // Clear buffer
	  s.resize(0);
	}
      else
	{
	  // Lookup locus name
	  ilocus = mlocus.find(t);
	  if (ilocus != mlocus.end())
	    s.push_back(ilocus->second);	  	  
	}

    }
  
  if (firsttime)
    printLOG(int2str(snpset.size()) 
	     + " sets read from [ " 
	     + par::setfile + " ] \n");
  
}



bool Plink::readClusterFile(bool verbose)
{

  checkFileExists(par::include_cluster_filename);
  ifstream CLST(par::include_cluster_filename.c_str(), ios::in);
  
  // Make map of family/individual IDs. Originally, set to not be
  // permuted: i.e. if an individual does not appear in the --within
  // file, then they will not be permuted. Otherwise, all individuals
  // are permuted. 

  map<string,Individual*> uid;
  map<string,Individual*>::iterator ii;
  map<string,int>::iterator ik;

  for (int i=0; i<sample.size(); i++)
    {
      uid.insert(make_pair(sample[i]->fid+"_"+sample[i]->iid,sample[i]));
      sample[i]->sol = -1;
    }
  
  int cnt=0;
  int k=0;
  kname.resize(0);
  kmap.clear();

  while (!CLST.eof())
    {
      string pfid, piid;
      string cluster;
      
      char cline[par::MAX_LINE_LENGTH];
      CLST.getline(cline,par::MAX_LINE_LENGTH,'\n');
      
      // convert to string
      string sline = cline;
      if (sline=="") continue;
      
      string buf; 
      stringstream ss(sline); 
      vector<string> tokens; 
      while (ss >> buf)
	tokens.push_back(buf);
      
      // Trying to read past last column 
      if (tokens.size() < 2+par::mult_clst) 
	{
	  if (! par::bmatch ) 
	    {
	      for (int i0=0; i0<tokens.size(); i0++)
		printLOG(tokens[i0]+" ");
	      printLOG("\n");
	    }
	  
	  CLST.close();
	  return false;
	} 
      
      pfid = tokens[0];
      piid = tokens[1];
      cluster = tokens[1 + par::mult_clst];

      // Find individual

      ii = uid.find(pfid+"_"+piid);
      
      if (ii != uid.end() )
	{
	  // Is this cluster on list?
	  ik = kmap.find(cluster);
	  if ( ik != kmap.end() ) // already on list
	    (ii->second)->sol = ik->second;
	  else // add to list
	    {	      
	      kmap.insert(make_pair(cluster,k));
	      (ii->second)->sol = k;
	      kname.push_back(cluster);
	      k++;
	    }
	  cnt++;
	}
    }

  CLST.close();

  // Assign to Family-like groups of pointers (e.g. to enable DFAM routine)

  klist.clear();
  for (int j=0;j<k;j++)
    klist.push_back( new Cluster );

  for (int i=0; i<n; i++)
    {
      if ( sample[i]->sol > -1 )
	klist[sample[i]->sol]->person.push_back(sample[i]);
    }
  
  if ( verbose )
    printLOG(int2str( cnt ) + " of " + int2str(n) +  
	     " individuals assigned to " 
	     + int2str(k) +  " cluster(s)\n");
  
  nk = k;
  
  return true;
}


bool Plink::readPhenoFile()
{

  // Assume binary trait unless we find out otherwise
  par::qt = false;
  par::bt = true; 
  
  printLOG( "Reading alternate phenotype from [ " 
	    + par::pheno_filename + " ] \n");

  checkFileExists(par::pheno_filename);
  ifstream PHE(par::pheno_filename.c_str(), ios::in);

  // Make map of family/individual IDs
  // and initially set all individuals as missing
  map<string,Individual*> uid;
  map<string,Individual*>::iterator ii;
  for (int i=0; i<sample.size(); i++)
    {
      uid.insert(make_pair(sample[i]->fid+"_"+sample[i]->iid,sample[i]));
      sample[i]->phenotype = -9;
      sample[i]->missing = true;
    }

  
  // Read in phenotype labels, if --pheno-name has been used

  if ( par::name_pheno != "" )
    {
      string pfid, piid, ph;
      char cline[par::MAX_LINE_LENGTH];
      PHE.getline(cline,par::MAX_LINE_LENGTH,'\n');
      string sline = cline;
      if (sline!="") 
	{
	  string buf; 
	  stringstream ss(sline); 
	  vector<string> tokens; 
	  while (ss >> buf)
	    tokens.push_back(buf);
	  
	  if ( tokens[0] != "FID" )
	    error("First header field must be FID");
	  if ( tokens[1] != "IID" )
	    error("Second header field must be IID");
	  
	  par::mult_pheno = -1;
	  for ( int i=2; i<tokens.size(); i++)
	    if ( tokens[i] == par::name_pheno )
	      par::mult_pheno = i - 1;
	  
	  if ( par::mult_pheno == -1 ) 
	    error("Did not find phenotype [ " 
		  + par::name_pheno 
		  + " ] in [ " 
		  + par::pheno_filename 
		  + " ]");
	}
    }

  int ccount = -1;

  while (!PHE.eof())
    {

      string pfid, piid, ph;
      char cline[par::MAX_LINE_LENGTH];
      PHE.getline(cline,par::MAX_LINE_LENGTH,'\n');
      
      // convert to string
      string sline = cline;
      if (sline=="") continue;

      string buf; 
      stringstream ss(sline); 
      vector<string> tokens; 
      while (ss >> buf)
	tokens.push_back(buf);

      if ( ccount < 0 ) 
	ccount = tokens.size();
      else if ( ccount != tokens.size() )
	error("Wrong number of columns in file [ "+par::pheno_filename
	      +" ] line :\n"+sline);
      
      if (tokens.size() < 2+par::mult_pheno) 
	{
	  if ( par::all_pheno )
	    {
	      printLOG("Processed all phenotypes\n");
	      return false;
	    }
	  else
	    error("Problem with [ "+par::pheno_filename
		  +" ] -- not enough columns :\n"+sline);
	}
      
      pfid = tokens[0];
      piid = tokens[1];

      // Skip header, if --pheno-name not used but should have been
      // (shouldn't matter in any case, as no individual with FID 
      // of "FID" should be allowed

      if ( pfid == "FID" || 
	   piid == "IID" )
	{
	  phenotype_name = tokens[1 + par::mult_pheno];
	  continue;
	}

      ph = tokens[1 + par::mult_pheno];
    
      // Are we using 0/1 coding?
      if (par::coding01) 
	{
	  if ( ph == "1" ) 
	    ph = "2";      
	  else if ( ph == "0" )
	    ph = "1";
	  else 
	    ph = "0";
	}
      
      ii = uid.find(pfid+"_"+piid);
      if (ii != uid.end() )
	{

	  (ii->second)->missing = true;
		      
	  // Only connect for non-missing phenotypes
	  if (ph != par::missing_phenotype)
	    {
	      (ii->second)->missing = false;
	    }

	  // Convert to double, checking for illegal values
	  if ( ! from_string<double>( (ii->second)->phenotype , ph , std::dec))
	    (ii->second)->missing = true;
	  
	  if (ph != par::missing_phenotype && 
	      ph != "0" && 
	      ph != "1" && 
	      ph != "2" ) 
	    {
	      par::qt = true;
	      par::bt = false; 	      
	    }
	  
	}
    }
  PHE.close();

  // If a binary trait, now make 0 missing also
  // i.e. if we never saw other than missing, 0, 1 or 2
  
  if (par::bt)
    for (int i=0; i<sample.size(); i++)
      if ( sample[i]->phenotype == 0 )
	sample[i]->missing = true;


  int new_nmissing=0;
  for (int i=0; i<sample.size(); i++)
    if ( !sample[i]->missing ) new_nmissing++;
  printLOG(int2str(new_nmissing) 
	   + " individuals with non-missing alternate phenotype\n");
  
  if (par::bt) 
    {

      if (par::coding01) 
	printLOG("Assuming a disease phenotype (0=unaff, 1=aff, other=miss)\n");
      else
	{
	  printLOG("Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)\n");
	  if (par::missing_phenotype!="0")
	    printLOG("Missing phenotype value is also " 
		     + par::missing_phenotype + "\n");
	}

      int ncase = 0;
      int ncontrol = 0;
      int nmissing = 0;
      for (int i=0; i<sample.size(); i++)
	{
	  if ( sample[i]->missing )
	    nmissing++;
	  else if ( sample[i]->phenotype == 1 )
	    ncontrol++;
	  else if ( sample[i]->phenotype == 2 )
	    ncase++;
	}
      printLOG(int2str(ncase)+" cases, "
	       +int2str(ncontrol)+" controls and "
	       +int2str(nmissing)+" missing\n");
    }
  else 
    {
      printLOG("Assuming a quantitative trait\n");
      printLOG("Missing phenotype value is " 
	       + par::missing_phenotype + "\n");
    }

  return true;    
}



void Plink::makePhenotype()
{

  // Implies a binary trait 
  par::qt = false;
  par::bt = true; 

  printLOG("Constructing a binary phenotype from [ "
	   +par::make_pheno_filename+" ]\n");

  if ( par::make_pheno_present ) 
    {
      map<string,int> imap;

      printLOG("All individuals not present will be set as controls\n");
      int cnt = 0;
      for (int i=0;i<n;i++)
	{
	  imap.insert(make_pair( sample[i]->fid +"_" + sample[i]->iid , i ));
	  sample[i]->phenotype = 1;
	  sample[i]->missing = false;
	  sample[i]->aff = false;
	}

      ifstream I1;
      while ( ! I1.eof() )
	{
	  vector<string> tokens = tokenizeLine(I1);
	  if ( tokens.size() < 2 )
	    continue;
	  string fiid = tokens[0] + "_" + tokens[1];
	  map<string,int>::iterator it = imap.find(fiid);
	  if ( it != imap.end() )
	    {
	      sample[ it->second ]->phenotype = 2;
	      sample[ it->second ]->aff = true;	      
	      ++cnt;
	    }
	}
      I1.close();
      printLOG( int2str( cnt ) + " individuals set as cases\n");
    }
  else
    {

      printLOG("Test value is [ " + par::make_pheno_value 
	       + " ] and missing value is [ " 
	       + par::missing_phenotype + " ]\n");
      
      // Swap filename as temporary cluster file
      
      string tmp_covar_file = par::include_cluster_filename;
      int tmp_mult_covar = par::mult_clst;
      
      par::include_cluster_filename = par::make_pheno_filename;
      par::mult_clst = 1;
      
      if (!readClusterFile())
	error("Problem reading filter file [ " + par::make_pheno_filename + " ]\n");
      
      // Put back the original covariate specificiation  
      par::include_cluster_filename = tmp_covar_file;
      par::mult_clst = tmp_mult_covar;
      
      int setCase = 0;
      int setControl = 0;
      int setMissing = 0;
      int notFound = 0;
      
      map<string,int>::iterator k = kmap.find( par::make_pheno_value );
      int value = k != kmap.end() ? k->second : -1;
      
      map<string,int>::iterator km = kmap.find( par::missing_phenotype );
      int missingValue = km != kmap.end() ? km->second : -1;
      
      for (int i=0; i<n; i++)
	sample[i]->missing = true;
      
      for (int i=0; i<n; i++)
	{
	  
	  Individual * person = sample[i];
	  
	  if ( person->sol == -1 ) 
	    {
	      person->missing = true;
	      person->phenotype = 0;
	      ++notFound;
	    }
	  else if ( person->sol == missingValue )
	    {
	      person->missing = true;
	      person->phenotype = 0;
	      ++setMissing;
	    }
	  else if ( person->sol == value ) 
	    {
	      person->missing = false;;
	      person->phenotype = 2;
	      ++setCase;
	    }
	  else
	    {
	      person->missing = false;;
	      person->phenotype = 1;
	      ++setControl;
	    }
	}
  
      printLOG("Set "+int2str(setCase)+" cases and "+int2str(setControl)+" controls, ");
      printLOG(int2str(setMissing)+" missing, "+int2str(notFound)+" not found\n");


      //////////////////////////////
      // Clear cluster values now
      
      nk=1;
      kmap.clear();
      kname.clear();
      for (int i=0; i<sample.size(); i++)
	sample[i]->sol = 0;
    }

}





bool Plink::readCovariateFile()
{

  // This will set individuals as missing 
  // if they have a missing value for the 
  // covariate, or do not appear in the file
  
  checkFileExists(par::covar_filename);
  ifstream COV(par::covar_filename.c_str(), ios::in);
  
  map<string,Individual*> uid;
  map<string,Individual*>::iterator ii;
  set<Individual*> hasCovariate;

  for (int i=0; i<sample.size(); i++)
    {
      uid.insert(make_pair(sample[i]->fid+"_"+sample[i]->iid,sample[i]));      
    }
  
  int nvalid=0;
  while (!COV.eof())
    {
      string pfid, piid, cov;
      
      char cline[par::MAX_LINE_LENGTH];
      COV.getline(cline,par::MAX_LINE_LENGTH,'\n');
      
      // convert to string
      string sline = cline;
      if (sline=="") continue;

      string buf; 
      stringstream ss(sline); 
      vector<string> tokens; 
      while (ss >> buf)
	tokens.push_back(buf);
      
      // Trying to read past last column 
      if (tokens.size() < 2+par::mult_covar) 
      {

	if (! par::qmatch )
	  {
	    for (int i0=0; i0<tokens.size(); i0++)
	      printLOG(tokens[i0]+" ");
	    printLOG("\n");
	  }

	COV.close();
        return false;
      } 
      
      pfid = tokens[0];
      piid = tokens[1];
      cov = tokens[1 + par::mult_covar];

      
      ii = uid.find(pfid+"_"+piid);
      if (ii != uid.end() )
	{
	  Individual * person = ii->second;


	  // Set covariate value
	  bool badValue = false;
	  if ( ! from_string<double>( person->covar, cov, std::dec ) )
	    badValue = true;
	  
	  // Note that we've seen a covariate for this individual
	  hasCovariate.insert(person);
	  
	  // Was this missing? 
	  if (cov == par::missing_phenotype || badValue )
	    {
	      person->missing = true;
	    }
	  else
	    nvalid++;
	  
	  
	}
    }
  COV.close();


  // Set to missing any individuals for who we did not see the covariate

  vector<Individual*>::iterator person = sample.begin();
  while ( person != sample.end() )
    {
      if ( hasCovariate.find( *person ) == hasCovariate.end() )
	(*person)->missing = true;
      person++;
    }
  
  
  printLOG("Reading covariate from [ " + par::covar_filename + " ] with ");
  printLOG("nonmissing values for "+int2str(nvalid)+" individuals\n");
  
  return true;
}

 
bool Plink::readCovListFile()
{

  // This will set individuals as missing if they have a missing value
  // for the covariate, or do not appear in the file
  
  
  ifstream COV(par::clist_filename.c_str(), ios::in);
  
  
  map<string,Individual*> uid;
  map<string,Individual*>::iterator ii;
  set<Individual*> hasCovariate;

  for (int i=0; i<sample.size(); i++)
    {
      uid.insert(make_pair(sample[i]->fid+"_"+sample[i]->iid,sample[i]));
    }

  // If need (for later selection) keep explicit track of what is missing
  map<Individual*, vector<bool> > isMissing;
  map<Individual*,bool> originalPersonMissingStatus;

  par::clist_number = -1;

  int nvalid=0;
  while (!COV.eof())
    {

      string pfid, piid;
      vector_t clist;

      char cline[par::MAX_LINE_LENGTH];
      COV.getline(cline,par::MAX_LINE_LENGTH,'\n');
      
      // convert to string
      string sline = cline;
      if (sline=="") continue;
      
      string buf; 
      stringstream ss(sline); 
      vector<string> tokens; 
      while (ss >> buf)
	tokens.push_back(buf);
      
      if ( par::clist_number < 0 ) 
	{
	  par::clist_number = tokens.size() - 2;
	  // Assign default headers (can be overwritten)
	  clistname.resize(par::clist_number);
	  for (int c=0; c<par::clist_number; c++)
	    clistname[c] = "COV"+int2str(c+1);	  
	}
      else if (tokens.size() != par::clist_number + 2 ) 
	{
	  printLOG("Line:\n"+sline+"\n");
	  COV.close();
	  return false;
	} 
      
      pfid = tokens[0];
      piid = tokens[1];
      
      ii = uid.find(pfid+"_"+piid);
      if (ii != uid.end() )
	{
	  Individual * person = ii->second;
	  
	  // Track individual covariate missing status
	  person->clistMissing.resize(par::clist_number);
	  
	  // Store original missingness status for this person
	  originalPersonMissingStatus.insert(make_pair( person, person->missing ));

	  vector<bool> missing_status;
	  
	  // Were any missing/bad values?
	  bool okay = true;

	  // Add covariate values to clist
	  person->clist.clear();
	  for (int c=2; c<par::clist_number+2; c++)
	    {
	      double t = 0;
	      if ( ! from_string<double>( t, tokens[c], std::dec ) )
		okay = false;
	      person->clist.push_back( t );
	    }
	  
	  // Note that we've seen a covariate for this individual
	  hasCovariate.insert(person);
	  	  
	  for (int c=0; c<par::clist_number; c++)
	    {
	      if ( tokens[c+2] == par::missing_phenotype )
		{
		  okay = false;
		  missing_status.push_back(false);
		  person->clistMissing[c] = true;
		}
	      else
		{
		  missing_status.push_back(true);
		  person->clistMissing[c] = false;
		}
	    }
	  
	  if (!okay)
	    person->missing = true;
	  else
	    nvalid++;

	  // Record, if we will use this below
	  if ( par::clist_selection )
	    isMissing.insert(make_pair( person, missing_status ));
	  
	}
      else if ( pfid == "FID" && piid == "IID" )
	{
	  // This is a header row -- read in covariate names
	  for (int c=0; c<par::clist_number; c++)
	    clistname[c] = tokens[c+2];
	}
    }
  COV.close();


  
  
  // Set to missing any individuals for who we did not see the covariate
  // But also fill their covariate list with missing values

  vector<Individual*>::iterator person = sample.begin();
  vector<bool> dummy_missing_status(n,false);
  while ( person != sample.end() )
    {
      if ( hasCovariate.find( *person ) == hasCovariate.end() )
	{
	  (*person)->missing = true;
	  (*person)->clist.clear();
	  (*person)->clist.resize(par::clist_number, -9 );
	  (*person)->clistMissing.clear();
	  (*person)->clistMissing.resize(par::clist_number, true );
	  
	  if ( par::clist_selection )
	    isMissing.insert(make_pair( (*person), dummy_missing_status ));
	}
      person++;
    }
  

  printLOG("Reading " 
	   + int2str(par::clist_number) 
	   + " covariates from [ " 
	   + par::clist_filename 
	   + " ] with ");
  printLOG("nonmissing values for "
	   +int2str(nvalid)
	   +" individuals\n");
  

  /////////////////////////////////////////////////////////
  // Do we actually want to keep all these covariates?
  
  if ( par::clist_selection_number || par::clist_selection_name )
    {
      vector<int> covlist;
      if ( par::clist_selection_number )
	{
	  NList nl(par::clist_number);
	  covlist = nl.deparseNumberList(par::clist_selection_string);
	}
      else
	{
	  map<string,int> mapping;
	  for (int c=0; c<par::clist_number; c++)
	    mapping.insert(make_pair( clistname[c],c));
	  NList nl(par::clist_number);
	  covlist = nl.deparseStringList(par::clist_selection_string,&mapping);
	}
      
      int nvalid = 0;
      for (int i=0; i<n; i++)
	{
	  Individual * person = sample[i];

	  vector_t tmp = person->clist;
	  vector<bool> tmpMissing = person->clistMissing;

	  person->clist.clear();
	  person->clistMissing.clear();

	  // Reset per-person missing code
	  person->missing = originalPersonMissingStatus.find( person )->second;
	  
	  vector<bool> missing_status = isMissing.find( person )->second;

	  bool okay = true;
	  for (int c=0; c<covlist.size(); c++)
	    {
	      person->clist.push_back( tmp[ covlist[c] ] );
	      person->clistMissing.push_back( tmpMissing[ covlist[c] ] );

	      if ( ! missing_status[covlist[c]] )
		{ 
		  person->missing = true;
		  okay = false; 
		}
	    }
	  
	  if ( okay ) nvalid++;
	}

      // Reset sample-wide values (names, number)

      vector<string> tmp = clistname;
      clistname.clear();
      for (int c=0; c<covlist.size(); c++)
	clistname.push_back( tmp[ covlist[c] ] );

      printLOG("Selected subset of " + int2str(covlist.size()) + " from " 
	       + int2str(par::clist_number) + " covariates\n");
      
      par::clist_number = covlist.size();
      
      printLOG("For these, nonmissing covariate values for "
	       +int2str(nvalid)+" individuals\n");
  

    }


  return true;
}


 
bool Plink::readMultiplePhenoFile()
{
  
  // This will set individuals as missing if they have a missing value
  // for the covariate, or do not appear in the file
  
  checkFileExists(par::multiple_phenotype_file);

  ifstream PHEFILE(par::covar_filename.c_str(), ios::in);
  
  map<string,Individual*> uid;
  map<string,Individual*>::iterator ii;
  set<Individual*> hasPhenotype;

  for (int i=0; i<sample.size(); i++)
    uid.insert(make_pair(sample[i]->fid+"_"+sample[i]->iid,sample[i]));
  
  // If need (for later selection) keep explicit track of what is missing

  map<Individual*, vector<bool> > isMissing;
  map<Individual*,bool> originalPersonMissingStatus;
  
  par::plist_number = -1;

  int nvalid=0;
  while (!PHEFILE.eof())
    {

      string pfid, piid;
      vector_t clist;

      vector<string> tokens = tokenizeLine( PHEFILE );

      if ( par::plist_number < 0 ) 
	{
	  par::plist_number = tokens.size() - 2;
	  // Assign default headers (can be overwritten)
	  plistname.resize(par::plist_number);
	  for (int c=0; c<par::plist_number; c++)
	    plistname[c] = "PHE"+int2str(c+1);	  
	}
      else if (tokens.size() != par::plist_number + 2 ) 
	{
	  printLOG("Line:\n");
	  for (int i=0;i<tokens.size(); i++) 
	    printLOG(tokens[i]+" ");
	  printLOG("\n");
	  PHEFILE.close();
	  return false;
	} 
      
      pfid = tokens[0];
      piid = tokens[1];
      
      ii = uid.find(pfid+"_"+piid);

      if (ii != uid.end() )
	{
	  Individual * person = ii->second;
	  
	  // Track individual phenotype missing status
	  person->plistMissing.resize(par::plist_number);
	  
	  // Store original missingness status for this person
	  originalPersonMissingStatus.insert(make_pair( person, person->missing ));

	  vector<bool> missing_status;
	  
	  // Were any missing/bad values?
	  bool okay = true;

	  // Add phenotype values to plist
	  person->plist.clear();
	  for (int c=2; c<par::plist_number+2; c++)
	    {
	      double t = 0;
	      if ( ! from_string<double>( t, tokens[c], std::dec ) )
		okay = false;
	      person->plist.push_back( t );
	    }
	  
	  // Note that we've seen a covariate for this individual
	  hasPhenotype.insert(person);
	  	  
	  for (int c=0; c<par::plist_number; c++)
	    {
	      if ( tokens[c+2] == par::missing_phenotype )
		{
		  okay = false;
		  missing_status.push_back(false);
		  person->plistMissing[c] = true;
		}
	      else
		{
		  missing_status.push_back(true);
		  person->plistMissing[c] = false;
		}
	    }
	  
	  if (!okay)
	    person->missing = true;
	  else
	    nvalid++;

	  // Record, if we will use this below
	  if ( par::plist_selection )
	    isMissing.insert(make_pair( person, missing_status ));
	  
	}
      else if ( pfid == "FID" && piid == "IID" )
	{
	  // This is a header row -- read in covariate names
	  for (int c=0; c<par::plist_number; c++)
	    plistname[c] = tokens[c+2];
	}
    }
  PHEFILE.close();
  
  // Set to missing any individuals for who we did not see the covariate
  // But also fill their covariate list with missing values
  
  vector<Individual*>::iterator person = sample.begin();
  vector<bool> dummy_missing_status(n,false);
  while ( person != sample.end() )
    {
      if ( hasPhenotype.find( *person ) == hasPhenotype.end() )
	{
	  (*person)->missing = true;
	  (*person)->plist.clear();
	  (*person)->plist.resize(par::plist_number, -9 );
	  (*person)->plistMissing.clear();
	  (*person)->plistMissing.resize(par::plist_number, true );
	  
	  if ( par::plist_selection )
	    isMissing.insert(make_pair( (*person), dummy_missing_status ));
	}
      person++;
    }
  

  printLOG("Reading " 
	   + int2str(par::plist_number) 
	   + " phenotypes from [ " 
	   + par::multiple_phenotype_file
	   + " ] with ");
  printLOG("nonmissing values for "
	   +int2str(nvalid)
	   +" individuals\n");
  

  /////////////////////////////////////////////////////////
  // Do we actually want to keep all these covariates?
  
  if ( par::plist_selection_number || par::plist_selection_name )
    {
      vector<int> phelist;
      if ( par::plist_selection_number )
	{
	  NList nl(par::plist_number);
	  phelist = nl.deparseNumberList(par::plist_selection_string);
	}
      else
	{
	  map<string,int> mapping;
	  for (int c=0; c<par::plist_number; c++)
	    mapping.insert(make_pair( plistname[c],c));
	  NList nl(par::plist_number);
	  phelist = nl.deparseStringList(par::plist_selection_string,&mapping);
	}
      
      int nvalid = 0;
      for (int i=0; i<n; i++)
	{
	  Individual * person = sample[i];
	  
	  vector_t tmp = person->plist;
	  vector<bool> tmpMissing = person->plistMissing;
	  
	  person->plist.clear();
	  person->plistMissing.clear();
	  
	  // Reset per-person missing code
	  person->missing = originalPersonMissingStatus.find( person )->second;
	  
	  vector<bool> missing_status = isMissing.find( person )->second;

	  bool okay = true;
	  for (int c=0; c<phelist.size(); c++)
	    {
	      person->plist.push_back( tmp[ phelist[c] ] );
	      person->plistMissing.push_back( tmpMissing[ phelist[c] ] );

	      if ( ! missing_status[phelist[c]] )
		{ 
		  person->missing = true;
		  okay = false; 
		}
	    }
	  
	  if ( okay ) nvalid++;
	}

      // Reset sample-wide values (names, number)
      
      vector<string> tmp = plistname;
      plistname.clear();
      for (int c=0; c<phelist.size(); c++)
	plistname.push_back( tmp[ phelist[c] ] );
      
      printLOG("Selected subset of " + int2str(phelist.size()) + " from " 
	       + int2str(par::plist_number) + " phenotypes\n");
      
      par::plist_number = phelist.size();
      
      printLOG("For these, nonmissing phenotype values for "
	       +int2str(nvalid)+" individuals\n");
  
    }


  return true;
}




void Plink::readSegmentFile(ifstream & SEG)
{
  
  // No need to skip header line(s) as no individuals should be called
  // "FID". Because of this, we can just concatenate multiple .segment
  // files, and not worry about repeating the headers

  map<string,Individual*> uid;
  for (int i=0; i<sample.size(); i++)
    uid.insert(make_pair(sample[i]->fid+"_"+sample[i]->iid,sample[i]));

  map<string,int> mlocus;
  for (int l=0; l<nl_all; l++)
    mlocus.insert(make_pair(locus[l]->name,l));

  map<string,Individual*>::iterator ii1;
  map<string,Individual*>::iterator ii2;

  map<string,int>::iterator il1;
  map<string,int>::iterator il2;
  
  int nseg=0;
  
  while (!SEG.eof())
    {
      string pfid1, piid1, pfid2, piid2;
      string snp1, snp2;
      string dummy;
	
      Segment s;
      
      // Contains pointers to individuals (p1, p2)
      // and ints for start/stop (SNP coding)      

      SEG >> pfid1 
	  >> piid1
	  >> pfid2
	  >> piid2 
	  >> dummy  // phenotype
	  >> dummy  // CHR
	  >> dummy  // BP1
	  >> dummy  // BP2
	  >> snp1
	  >> snp2
	  >> dummy  // NSNP
	  >> dummy; // KB


      bool okay = true;

      // Attached individuals
      ii1 = uid.find(pfid1+"_"+piid1);
      ii2 = uid.find(pfid2+"_"+piid2);
      if ( ii1 != uid.end() && ii2 != uid.end() )
	{
	  s.p1 = ii1->second;
	  s.p2 = ii2->second;
	}
      else okay = false;
	
      // Attached bounding SNPs
      il1 = mlocus.find(snp1);
      il2 = mlocus.find(snp2);
      if ( il1 != mlocus.end() && il2 != mlocus.end() )
	{
	  s.start = il1->second;
	  s.finish = il2->second;
	}
      else okay = false;
      
      if ( okay ) 
	  if ( locus[s.finish]->bp - locus[s.start]->bp 
	       < par::segment_length || 
	       s.finish - s.start + 1 < par::segment_snp ) okay = false;
      
      // Add to list
      if (okay)
      {
	  segment.push_back(s);
	  nseg++;
	}

    }
  
  printLOG("Read " + int2str(nseg) + " valid segments\n");

}



void Plink::readSegmentFileMinimal(ifstream & SEG)
{
  // Format: (no header)
  
  // {n1, n2} , {s1,s2}, {s1,s2} , {s1,s2} , {-1,-1} \n 
  
  // n1 = int for individual 1 position in sample
  // n2 = int for individual 2 position in sample
  // s1 = segment start (first SNP int position)
  // s2 = segment finish

  // *** Health warning : the sample must be identical, or else markers
  // and individuals will be out... The verbose format avoids this problem, 
  // but will generate much larger files ***

  int nseg=0;
  
  while (!SEG.eof())
    {
      int p1, p2;

      Segment s;
      
      // Contains pointers to individuals (p1, p2)
      // and ints for start/stop (SNP coding)      

      // For this pair;

      SEG >> p1 >> p2;
      
      if ( p1 < 0 || p2 > n )
	error("Problem with segment file (minimal format)\n");

      s.p1 = sample[p1];
      s.p2 = sample[p2];
      
      // Read segments
      while (true)
	{
	  
	  SEG >> s.start >> s.finish;
	  
	  // Move to next pair? 	  
	  
	  if ( s.start == -1 ) 
	    break;
	  
	  // Segment long enough?
	  
	  if ( locus[s.finish]->bp - locus[s.start]->bp 
	       >= par::segment_length &&
	       s.finish - s.start + 1 >= par::segment_snp ) 
	    {
	      segment.push_back(s);
	      nseg++;
	    }
	  
	  // Read next segment for this pair
	}
      
      // Read next pair, while not EOF
    }
    
  printLOG("Read " + int2str(nseg) + " valid segments\n");

}


void Plink::readConditioningList()
{

  checkFileExists(par::conditioning_snps_file);
  ifstream COV(par::conditioning_snps_file.c_str(), ios::in);

  printLOG("Reading list of conditioning SNPs from [ " 
	   + par::conditioning_snps_file + " ]\n");

  int c=0;
  int c2=0;
  while (!COV.eof())
    {
      string snp;
      COV >> snp;
      if (snp=="") 
	break;
      int x = getMarkerNumber( *this , snp );
      if (x>=0)
	{
	  conditioner.push_back( getMarkerNumber( *this , snp ) );
	  conditioner_mask.push_back( false );
	  c++;
	}
      c2++;
    }

  printLOG("Using " + int2str(c) + " of "
	   +int2str(c2)
	   +" specified conditioning SNPs\n");
  COV.close();

}


void Plink::readMapFile(string filename,
			vector<bool> & include,
			vector<int> & include_pos,
			int & nl_actual )
{
  
  // chromosome code
  // SNP identifier   
  // cM / M      
  // Base Position (-ve implies exclude)
  
  vector<Locus> ordered;
  
  ifstream MAP;
  MAP.open( filename.c_str());
  MAP.clear();

  int c=0;
  while(!MAP.eof())
    {
           
      string chr;
      long int inc;

      char cline[256];
      MAP.getline(cline,256,'\n');

      // convert to string
      string sline = cline;
      if (sline=="") 
	continue;
      
      string buf; 
      stringstream ss(sline); 
      vector<string> tokens; 
      while (ss >> buf)
	tokens.push_back(buf);
      
      if (tokens.size() == 0)
	continue;
     
      // Is this line a comment?      
      if (tokens[0].substr(0,1)=="#")
	{
	  continue;
	}
      
      if ( par::map3 && tokens.size() != 3 ) 
	error("Problem with MAP file line:\n"+sline);
      else if ( (!par::map3) && tokens.size() != 4 )
	error("Problem with MAP file line:\n"+sline);
      
      
      Locus * loc;

      if ( ! par::load_gvar )
	loc = new Locus;
      else
	{
	  Variant * vloc = new Variant;
	  loc = (Locus*)vloc; 
	}

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

  MAP.clear();
  MAP.close();
   
  printLOG(int2str(locus.size()) + " (of " + int2str(include.size()) + 
	   ") markers to be included from [ " + filename + " ]\n");
  
  if ( par::load_gvar ) 
    printLOG("Read in as generic variants, rather than SNPs\n");
  
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

  nl_actual = locus.size();
  
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
  
  if (par::SNP_major && ! par::load_gvar )
    {
      for (int i=0; i<nl_actual; i++)
	{
	  CSNP * newlocus = new CSNP;
	  SNP.push_back(newlocus);
	}
    }
  
}


void Plink::readFamFile(string filename)
{

  //////////////////////////////
  // Read pedigree information

  // Initially, assume a binary trait
  par::qt = false;
  par::bt = true; 

  printLOG("Reading pedigree information from [ " 
	   + filename + " ] \n");

  checkFileExists( filename );

  ifstream PED;
  PED.open(filename.c_str());
  PED.clear();

  vector<Individual*> ambiguous;
  
  int nmale = 0;
  int nfemale = 0;
  int nambig = 0;

  int c=0;
  while(!PED.eof())
    {

      Individual * person = new Individual;

      // First 6 obligatory fields
      string phenotype;

      PED >> person->fid 
	  >> person->iid 
	  >> person->pat
	  >> person->mat
	  >> person->sexcode
	  >> phenotype;
      
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
      if (person->fid=="") 
	{
	  delete person;
	  break;
	}

      // Check for reserved family ID code
      if ( person->fid=="FID" )
	error("FID is a reserved ID... please select a different family ID");

      // Check sex
      if (person->sexcode=="1")
	{
	  person->sex = true; // male
	  nmale++;
	}
      else if (person->sexcode=="2")
	{
	  person->sex = false;  // female (default)
	  nfemale++;
	}
      else 
	{
	  ambiguous.push_back(person);
	  nambig++;
	  if (!par::ignore_missing_sex)
	    person->missing = true;
	}

      
      ///////////////
      // A non-founder?
      
      person->founder = (person->pat == "0" && person->mat == "0") ? true : false;
      
      //////////////////////////////
      // Test for quantitative trait
      
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
      
      
      // Increase person counter
      c++;
      
      // Add individual to list
      sample.push_back(person);
      
    }
    
  PED.clear();
  PED.close();
  

  // If a binary trait, now make 0 missing also
  // i.e. if we never saw other than missing, 0, 1 or 2
  
  if (par::bt)
    for (int i=0; i<sample.size(); i++)
      if ( sample[i]->phenotype == 0 )
	sample[i]->missing = true;

  
  printLOG(int2str(c)+" individuals read from [ " 
	   + filename + " ] \n");
  int nm=0;
  for (int i=0;i<sample.size();i++)
    if(!sample[i]->missing) nm++;
  printLOG(int2str(nm) 
	   + " individuals with nonmissing phenotypes\n");
  
  if (par::bt) 
    {

      if (par::coding01) 
	printLOG("Assuming a disease phenotype (0=unaff, 1=aff, other=miss)\n");
      else
	{
	  printLOG("Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)\n");
	  if (par::missing_phenotype!="0")
	    printLOG("Missing phenotype value is also " 
		     + par::missing_phenotype + "\n");
	}

      int ncase = 0;
      int ncontrol = 0;
      int nmissing = 0;

      for (int i=0; i<sample.size(); i++)
	if ( sample[i]->missing )
	  nmissing++;
	else if ( sample[i]->phenotype == 1 )
	  ncontrol++;
	else if ( sample[i]->phenotype == 2 )
	  ncase++;
      printLOG(int2str(ncase)+" cases, "
	       +int2str(ncontrol)+" controls and "
	       +int2str(nmissing)+" missing\n");
      
    }
  else 
    {
      printLOG("Assuming a quantitative trait\n");
      printLOG("Missing phenotype value is " + par::missing_phenotype + "\n");
    }

  // Display sex counts
  printLOG(int2str(nmale)+" males, "+int2str(nfemale)
	   +" females, and "+int2str(nambig)+" of unspecified sex\n");


  // Display list of ambiguously-sexed individuals?
  if (ambiguous.size()>0)
    {
      printLOG("Warning, found "+int2str(ambiguous.size())
	       +" individuals with ambiguous sex codes\n");
      if (!par::ignore_missing_sex)
	printLOG("These individuals will be set to missing ( or use --allow-no-sex )\n");      
      string f = par::output_file_name + ".nosex";
      printLOG("Writing list of these individuals to [ "+f+" ]\n");
      ofstream AMB;
      AMB.open(f.c_str(), ifstream::out);
      for (int i=0; i<ambiguous.size(); i++)
	AMB << ambiguous[i]->fid << "\t" << ambiguous[i]->iid << "\n";
      AMB.close();      
      ambiguous.clear();
    }


}


void Plink::readHomozygSegmentFile(ifstream & SEG)
{
  
  // No need to skip header line(s) as no individuals should be called
  // "FID". Because of this, we can just concatenate multiple .segment
  // files, and not worry about repeating the headers
  
  map<string,Individual*> uid;
  for (int i=0; i<sample.size(); i++)
    uid.insert(make_pair(sample[i]->fid+"_"+sample[i]->iid,sample[i]));

  map<string,int> mlocus;
  for (int l=0; l<nl_all; l++)
    mlocus.insert(make_pair(locus[l]->name,l));

  map<string,Individual*>::iterator ii;
  map<string,int>::iterator il1;
  map<string,int>::iterator il2;
  
  int nseg=0;
  
  while (!SEG.eof())
    {
      string pfid, piid;
      string snp1, snp2;
      string dummy;
      
      Segment s;
      
      // Contains pointers to individuals (p1, p2)
      // and ints for start/stop (SNP coding)      

      SEG >> pfid 
	  >> piid
	  >> dummy   // phenotype
	  >> dummy   // CHR
	  >> snp1
	  >> snp2
	  >> dummy   // BP1
	  >> dummy   // BP2
	  >> dummy   // NSNP
	  >> dummy   // KB
	  >> dummy   // DENSITY
	  >> dummy   // PHOM
	  >> dummy;  // PHET

      bool okay = true;

      // Attached individuals
      ii = uid.find(pfid+"_"+piid);
      if ( ii != uid.end() ) 
	s.p1 = s.p2 = ii->second;
      else 
	okay = false;
	
      // Attached bounding SNPs
      il1 = mlocus.find(snp1);
      il2 = mlocus.find(snp2);
      if ( il1 != mlocus.end() && il2 != mlocus.end() )
	{
	  s.start = il1->second;
	  s.finish = il2->second;
	}
      else okay = false;
      
      if ( okay ) 
      {
	  if ( locus[s.finish]->bp - locus[s.start]->bp 
	       < ( par::homo_run_length_kb * 1000) || 
	       s.finish - s.start + 1 < par::homo_run_length_snps ) 
	      okay = false;
      }

      // Add to list
      if (okay)
      {
	  segment.push_back(s);
	  nseg++;
      }
      
    }
  
  printLOG("Read " + int2str(nseg) + " valid segments\n");

}


void Plink::readStdIn()
{

  // Read a PED/MAP file from standard input
  

}


void Plink::updateMapFile()
{
  
  map<string,int> mlocus;
  for (int l=0; l<nl_all; l++)
    mlocus.insert(make_pair( locus[l]->name, l ));

  if ( par::update_cm )
    printLOG("Reading new cM/M positions from [ "
	     + par::update_mapfile + " ]\n");
  else if ( par::update_chr )
    printLOG("Reading new chromosome positions from [ "
	     + par::update_mapfile + " ]\n");
  else if ( par::update_name )
    printLOG("Reading new SNP labels from [ "
	     + par::update_mapfile + " ]\n");
  else
    printLOG("Reading new physical positions from [ "
	     + par::update_mapfile + " ]\n");
  
  checkFileExists( par::update_mapfile );

  ifstream MAPIN;
  
  MAPIN.open( par::update_mapfile.c_str(), ios::in );

  int num_found = 0;
  int num_notfound = 0;

  set<string> done;
  set<string> names;
  bool nameWarning = false;

  while ( ! MAPIN.eof() )
    {
      
      string snp;
      double value;
      string svalue;
      
      if ( par::update_chr || par::update_name )
	MAPIN >> snp >> svalue;
      else
	MAPIN >> snp >> value;
      
      if ( snp == "" )
	continue;
      
      map<string,int>::iterator il = mlocus.find( snp );
      
      set<string>::iterator sl = done.find( snp );      

      if ( sl != done.end() )
	{
	  error(snp+" seen more than once in [ "+par::update_mapfile+" ]\n");
	}
      
      
      if ( il != mlocus.end() )
	{

	  done.insert( snp );

	  if ( par::update_name ) 
	    {
	      if ( names.find( svalue ) != names.end() )
		nameWarning = true;
	      names.insert( svalue );
	    }
	  
	  ++num_found;
	  
	  if ( par::update_cm ) 
	    locus[ il->second ]->pos = value;
	  else if ( par::update_chr )
	    locus[ il->second ]->chr = getChromosomeCode(svalue);
	  else if ( par::update_name )
	    locus[ il->second ]->name = svalue;
	  else
	    locus[ il->second ]->bp = (int)value;	  
	}
      else
	++num_notfound;
      
    }

  printLOG( int2str( num_found ) + " SNP positions read and updated\n");
  printLOG( int2str( nl_all - done.size() ) + " in data but not in [ " + par::update_mapfile + " ]\n");
  if ( num_notfound > 0 )
    printLOG(int2str(num_notfound)
	     + " in [ " + par::update_mapfile + " ] but not in data\n");
  

  if ( par::update_name && nameWarning ) 
    printLOG("Warning -- duplicated SNP names found in update\n");
 
  // Check if file needs re-ordering
  
  bool allInOrder = true;
  
  for (int l=1; l<nl_all; l++)
    {
      if ( locus[l]->chr == locus[l-1]->chr ) 
	{
	  if ( par::update_cm ) 
	    {
	      if ( locus[l]->pos < locus[l-1]->pos )
		allInOrder = false;
	    }
	  else
	    {
	      if ( locus[l]->bp < locus[l-1]->bp )
		allInOrder = false;
	    }
	}
      else if ( locus[l]->chr < locus[l-1]->chr ) 
	allInOrder = false;

      if ( ! allInOrder ) 
	break;
    }

  if ( ! allInOrder )
    printLOG("*** Implicit order changed from re-mapping ***\n");
  
}



void Plink::updateAlleles()
{
  
  map<string,int> mlocus;
  for (int l=0; l<nl_all; l++)
    mlocus.insert(make_pair( locus[l]->name, l ));

  printLOG("Reading new allele codes from [ " + par::update_allele_file + " ]\n");
  
  checkFileExists( par::update_allele_file );

  ifstream MAPIN;
  
  MAPIN.open( par::update_allele_file.c_str(), ios::in );

  int num_found = 0;
  int num_notfound = 0;
  int num_prob = 0;

  set<string> done;
  set<int> probs;

  while ( ! MAPIN.eof() )
    {
      
      string snp;
      string old1, old2, new1, new2;
      
      MAPIN >> snp >> old1 >> old2 
	    >> new1 >> new2;
      
      if ( snp == "" )
	continue;
      
      map<string,int>::iterator il = mlocus.find( snp );
      
      set<string>::iterator sl = done.find( snp );      
      
      if ( sl != done.end() )
	{
	  error(snp + " seen more than once in [ "+par::update_allele_file+" ]\n");
	}
      
      
      bool success = false;
      
      if ( il != mlocus.end() )
	{
	  done.insert( snp );
	  
	  ++num_found;	  
	  
	  if ( locus[ il->second ]->allele1 == old1 || 
	       locus[ il->second ]->allele2 == old2 )
	    {
	      locus[ il->second ]->allele1 = new1;
	      locus[ il->second ]->allele2 = new2;
	    }
	  else if ( locus[ il->second ]->allele1 == old2 || 
		    locus[ il->second ]->allele2 == old1 )
	    {
	      locus[ il->second ]->allele1 = new2;
	      locus[ il->second ]->allele2 = new1;
	    }
	  else  
	    {
	      ++num_prob;
	      probs.insert( il->second );
	    }
	}
      else
	++num_notfound;
      
    }

  printLOG( int2str( num_found-num_prob ) + " SNPs found and allele codes updated\n");
  if( nl_all - done.size() > 0 )
    printLOG( int2str( nl_all - done.size() ) + " in data but not in [ " + par::update_allele_file + " ]\n");

  if ( num_notfound > 0 )
    printLOG(int2str(num_notfound)
	     + " in [ " + par::update_allele_file + " ] but not in data\n");

  if ( num_prob > 0 )
    {
      printLOG(int2str(num_prob)
	       + " SNPs with allele conflicts listed in [ " 
	       + par::output_file_name + ".allele.no.snp ]\n");
      ofstream O( ( par::output_file_name + ".allele.no.snp" ).c_str() , ios::out);
      set<int>::iterator i = probs.begin();
      while ( i != probs.end() )
	{
	  O << locus[ *i ]->name << "\t"
	    << locus[ *i ]->allele1 << "\t"
	    << locus[ *i ]->allele2 << "\n";
	  
	  ++i;
	}
      O.close();
      
    }
  
}


void Plink::updateFamFile()
{
  
  map<string,int> mpeople;
  for (int i=0; i<n; i++)
    mpeople.insert(make_pair( sample[i]->fid+"_"+sample[i]->iid , i ));
  
  ifstream FAM_ID, FAM_PAR, FAM_SEX, FAM_PHE;

  if ( par::update_ids )
    {
      printLOG("Reading new FIDs and IIDs from [ "
	     + par::update_ids_file + " ]\n");
      checkFileExists( par::update_ids_file );
      FAM_ID.open( par::update_ids_file.c_str(), ios::in );
    }

  
  if ( par::update_ids )
    {
      int num_found = 0;
      int not_found = 0; 
      while (! FAM_ID.eof() )
	{
	  vector<string> tokens = tokenizeLine( FAM_ID );
	  if ( tokens.size() == 0 )
	    continue;
	  if ( tokens.size() != 4 ) 
	    error("Problem with line in --update-ids file: expects 4 columns per row");
	  
	  string index = tokens[0] + "_" + tokens[1];	  
	  if ( mpeople.find( index ) != mpeople.end() )
	    {
	      ++num_found;
	      Individual * p = sample[ mpeople.find(index)->second ];
	      p->fid = tokens[2];
	      p->iid = tokens[3];
	    }
	  else
	    ++not_found;
	}
      printLOG( int2str(num_found) + " individuals found, " 
		+ int2str(not_found) + " not in sample\n");
      FAM_ID.close();
    }
  
  if ( par::update_sex )
    {
      printLOG("Reading new sex codes from [ "
	       + par::update_sex_file + " ]\n");
      checkFileExists( par::update_sex_file );
      FAM_SEX.open( par::update_sex_file.c_str(), ios::in );
    }
  
  if ( par::update_sex )
    {
      int num_found = 0;
      int not_found = 0;
      while (! FAM_SEX.eof() )
	{
	  vector<string> tokens = tokenizeLine( FAM_SEX );
	  if ( tokens.size() == 0 )
	    continue;
	  if ( tokens.size() != 3 ) 
	    error("Problem with line in --update-sex file: expects 3 columns per row");
	  
	  string index = tokens[0] + "_" + tokens[1];	  
	  if ( mpeople.find( index ) != mpeople.end() )
	    {
	      ++num_found;
	      Individual * p = sample[ mpeople.find(index)->second ];
	      p->sexcode = tokens[2];
	      
	      if (p->sexcode=="1")
		p->sex = true; // male
	      else if (p->sexcode=="2")
		p->sex = false;  // female
	      else if (!par::ignore_missing_sex)
		p->missing = true;
	    }
	  else
	    ++not_found;
	}
      printLOG( int2str(num_found) + " individuals found, " 
		+ int2str(not_found) + " not in sample\n");
      FAM_SEX.close();
    }

    
    if ( par::update_parents )
    {
      printLOG("Reading new parental codes from [ "
	       + par::update_parents_file + " ]\n");
      checkFileExists( par::update_parents_file );
      FAM_PAR.open( par::update_parents_file.c_str(), ios::in );
    }
    
    if ( par::update_parents )
      {
	int num_found = 0;
	int not_found = 0;
	while (! FAM_PAR.eof() )
	{
	  vector<string> tokens = tokenizeLine( FAM_PAR );
	  if ( tokens.size() == 0 )
	    continue;
	  if ( tokens.size() != 4 ) 
	    error("Problem with line in --update-parents file: expects 4 columns per row");
	  
	  string index = tokens[0] + "_" + tokens[1];	  
	  if ( mpeople.find( index ) != mpeople.end() )
	    {
	      ++num_found;
	      Individual * p = sample[ mpeople.find(index)->second];
	      p->pat = tokens[2];
	      p->mat = tokens[3];
	    }
	  else
	    ++not_found;
	}
      printLOG( int2str(num_found) + " individuals found, " 
		+ int2str(not_found) + " not in sample\n");
      FAM_PAR.close();
    }

    if ( par::update_pheno )
    {
      printLOG("Reading phenotypes to update from [ "
	       + par::update_pheno_file + " ]\n");
      checkFileExists( par::update_pheno_file );
      FAM_PHE.open( par::update_pheno_file.c_str(), ios::in );
    }
    
    if ( par::update_pheno )
      {
	int num_found = 0;
	int not_found = 0;
	while (! FAM_PHE.eof() )
	{
	  vector<string> tokens = tokenizeLine( FAM_PHE );
	  if ( tokens.size() == 0 )
	    continue;
	  if ( tokens.size() != 3 ) 
	    error("Problem with line in --update-pheno file: expects 3 columns per row");
	  
	  string index = tokens[0] + "_" + tokens[1];	  
	  if ( mpeople.find( index ) != mpeople.end() )
	    {
	      ++num_found;
	 
	      Individual * p = sample[ mpeople.find(index)->second];
	      
	      if (par::coding01) 
		{
		  if ( tokens[2] == "1" ) 
		    tokens[2] = "2";      
		  else if ( tokens[2] == "0" )
		    tokens[2] = "1";
		  else 
		    tokens[2] = "0";
		}

	      if ( tokens[2] == par::missing_phenotype)
		p->missing = true;
	      else
		{
		  if (  ! from_string<double>( p->phenotype, tokens[2], std::dec ) )
		    p->missing = true;
		  else
		    {
		      if ( tokens[2] != "0" && 
			   tokens[2] != "1" && 
			   tokens[2] != "2" ) 
			{
			  par::qt = true;
			  par::bt = false; 
			}
		    }
		}
	    }
	  else
	    ++not_found;
	}

	printLOG( int2str(num_found) + " individuals found, " 
		  + int2str(not_found) + " not in sample\n");
	FAM_PHE.close();

	// Do we need to recode binary phenotypes?
	if (par::bt)
	  affCoding(*this);

      }
    
}

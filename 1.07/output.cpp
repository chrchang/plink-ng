

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
#include <climits>
#include <fstream>
#include <sstream>
#include <bitset>
#include <map>

#include "plink.h"
#include "options.h"
#include "perm.h"
#include "sets.h"
#include "helper.h"
#include "nlist.h"

extern ofstream LOG;

void Plink::printLOG(string s)
{
  LOG << s;
  LOG.flush();
  
  if (!par::silent)
    {
      cout << s;
      cout.flush();
    }
}

void Plink::display_indivReport()
{
  Individual * p1 = NULL;
  int i1;
  
  for (int i=0; i<n; i++)
    {
      if ( sample[i]->fid == par::indiv_report_fid && 
	   sample[i]->iid == par::indiv_report_iid )
	{ 
	  p1 = sample[i]; 
	  i1 = i; 
	  break;
	}
    }
  
  if ( p1 == NULL )
    error("Problem finding individual indicated in --report\n");
  
  printLOG("\nReport for individual [ " + p1->fid + " " + p1->iid + " ]\n\n");
  
  for ( int l = 0 ; l < nl_all ; l++)
    {
      stringstream s2;
      s2 << setw(4) << locus[l]->chr << " " 
	 << setw(par::pp_maxsnp) << locus[l]->name << " "
	 << setw(6) << genotype(*this,i1,l) << "\n";
      printLOG( s2.str() );
    }
  
}

void Plink::displayGenomePV()
{

  ofstream PLO;
  string f = par::output_file_name + ".plink2";
  PLO.open(f.c_str(),ios::out);
  
  printLOG("Writing genome-wide corrected PLINK results to [ " +f +" ] \n");
  
  for (int i=0; i<original.size();i++)
    {
      for (int j=0; j<original[i].size(); j++)
	{
	  double pv=0;
	  for (int k=0; k<maxr2.size(); k++)
	    if (maxr2[k] >= original[i][j]) pv++;
	  PLO << "G " 
	      << double(pv+1)/double(par::replicates+1)
	      << "\n";
	}
    }
  PLO.close();
}


void Plink::display_pairList()
{

  if (par::SNP_major)
    SNP2Ind();
  
  string f = par::output_file_name + ".plist";
  
  printLOG("Writing pair-list file to [ " + f + " ] \n");
  ofstream PED(f.c_str(), ios::out);
  PED.clear();
  PED.precision(4);

  // Find individuals
  Individual * p1 = NULL;
  Individual * p2 = NULL;
  int i1, i2;

  for (int i=0; i<n; i++)
    {
      if ( sample[i]->fid == par::plist_fid1 && sample[i]->iid == par::plist_iid1 )
	{ p1 = sample[i]; i1 = i; }
      if ( sample[i]->fid == par::plist_fid2 && sample[i]->iid == par::plist_iid2 )
	{ p2 = sample[i]; i2 = i; } 
    }
  
  if (p1 == NULL || p2 == NULL)
    error("Problem finding individuals indicated in --plist\n");


  PED << setw(4) << "CHR" << " " 
      << setw(par::pp_maxsnp) << "SNP" << " "
      << setw(14) << "BP" << " "
      << setw(par::pp_maxfid+par::pp_maxiid) << (p1->fid+"/"+p1->iid) << " "
      << setw(par::pp_maxfid+par::pp_maxiid) << (p2->fid+"/"+p2->iid) << " "
      << setw(4) << "A1" << " "
      << setw(8) << "MAF" << " "
      << setw(4) << "IBS" << "\n";
  
  for (int l=0; l<nl_all; l++)  
    {
      PED << setw(4) << locus[l]->chr << " "
	  << setw(par::pp_maxsnp) << locus[l]->name << " "
	  << setw(14) << locus[l]->bp << " "
	  << setw(par::pp_maxfid+par::pp_maxiid) << genotype(*this,i1,l) << " "
	  << setw(par::pp_maxfid+par::pp_maxiid) << genotype(*this,i2,l) << " "
	  << setw(4) << locus[l]->allele1 << " "
	  << setw(8) << locus[l]->freq << " ";
      
      bool a1 = p1->one[l];
      bool a2 = p1->two[l];

      bool b1 = p2->one[l];
      bool b2 = p2->two[l];
	
      bool ibs0 = false;
      bool ibs1 = false;
      bool miss = false;
      
      if ( a1 == a2 && 
	   b1 == b2 && 
	   a1 != b1  ) ibs0 = true;
      else if ( a1 && !(a2) ) miss = true;
      else if ( b1 && !(b2) ) miss = true;
      else if ( a1 != b1 ||
		a2 != b2 ) ibs1 = true;

      PED << setw(4);
      if (ibs0) PED << "0" << "\n";
      else if (ibs1) PED << "1" << "\n";
      else if (miss) PED << "NA" << "\n";
      else PED << "2" << "\n";
      
    }
  
  PED.close();

}

void Plink::display_listByAllele()
{

  if (!par::SNP_major)
    Ind2SNP();

  // Create an output file that lists one allele per line
  // N = number of individuals FID1 IID2 FID2 IID2 ... FIDN IIDN 
  
  // SNP ALLELE N { Individual list ( 2N entries) } 

  string f = par::output_file_name + ".list";
  
  printLOG("Writing recoded list file to [ " + f + " ] \n");
  ofstream PED(f.c_str(), ios::out);
  PED.clear();

  vector<CSNP*>::iterator s = SNP.begin();
  vector<Locus*>::iterator loc = locus.begin();

  while ( s != SNP.end() )
    {
      // Genotype 11
      PED << (*loc)->chr << par::recode_delimit
	  << (*loc)->name << par::recode_delimit
	  << (*loc)->allele1 << (*loc)->allele1 ;
      
      vector<bool>::iterator i1 = (*s)->one.begin();
      vector<bool>::iterator i2 = (*s)->two.begin();
      vector<Individual*>::iterator gperson = sample.begin();
      
      while ( gperson != sample.end() )
	{
	  if ( (!*i1) && (!*i2) ) 
	    PED << par::recode_delimit << (*gperson)->fid 
		<< par::recode_delimit << (*gperson)->iid ;
	  i1++;
	  i2++;
	  gperson++;
	}
      PED << "\n";
      
      // Genotype 12
      PED << (*loc)->chr << par::recode_delimit
	  << (*loc)->name << par::recode_delimit
	  << (*loc)->allele1 << (*loc)->allele2 ;
      i1 = (*s)->one.begin();
      i2 = (*s)->two.begin();
      gperson = sample.begin();
      while ( gperson != sample.end() )
	{
	  if ( (!*i1) && *i2 ) 
	    PED << par::recode_delimit << (*gperson)->fid 
		<< par::recode_delimit << (*gperson)->iid ;
	  i1++;
	  i2++;
	  gperson++;
	}
      PED << "\n";

      // Genotype 22
      PED << (*loc)->chr << par::recode_delimit
	  << (*loc)->name << par::recode_delimit
	  << (*loc)->allele2 << (*loc)->allele2 ;
      i1 = (*s)->one.begin();
      i2 = (*s)->two.begin();
      gperson = sample.begin();
      while ( gperson != sample.end() )
	{
	  if ( *i1 && *i2 ) 
	    PED << par::recode_delimit << (*gperson)->fid 
		<< par::recode_delimit << (*gperson)->iid ;
	  i1++;
	  i2++;
	  gperson++;
	}
      PED << "\n";

      // Genotype 00
      PED << (*loc)->chr << par::recode_delimit
	  << (*loc)->name << par::recode_delimit
	  << "00";
      i1 = (*s)->one.begin();
      i2 = (*s)->two.begin();
      gperson = sample.begin();
      while ( gperson != sample.end() )
	{
	  if ( *i1 && !*i2 ) 
	    PED << par::recode_delimit << (*gperson)->fid 
		<< par::recode_delimit << (*gperson)->iid ;
	  i1++;
	  i2++;
	  gperson++;
	}
      PED << "\n";
      
      
      // Next SNP
      s++;
      loc++;
    }
  
  PED.close();

}

void make2LTable(ofstream & TWOL, Plink & P, int m1, int m2, vector_t count, bool percent)
{
  vector_t margA(4);
  vector_t margB(4);
  double total = 0;
  
  for (int i=0; i<16; i++)
    total += count[i];

  if (percent)
    {
      for (int i=0; i<16; i++)
	count[i] /= total;
      total = 1;
    }

  margA[0] = count[0]+count[1]+count[2]+count[3];
  margA[1] = count[4]+count[5]+count[6]+count[7];
  margA[2] = count[8]+count[9]+count[10]+count[11];
  margA[3] = count[12]+count[13]+count[14]+count[15];

  margB[0] = count[0]+count[4]+count[8]+count[12];
  margB[1] = count[1]+count[5]+count[9]+count[13];
  margB[2] = count[2]+count[6]+count[10]+count[14];
  margB[3] = count[3]+count[7]+count[11]+count[15];


  TWOL << setw(par::pp_maxsnp) << " " << " "
       << setw(4) << " " << " "
       << "   " << par::twolocus_snp2 << "\n";
  
  TWOL << setw(par::pp_maxsnp) << " " << " "
       << setw(4) << " " << " "
       << setw(6) << (P.locus[m2]->allele1+"/"+P.locus[m2]->allele1) << " "
       << setw(6) << (P.locus[m2]->allele1+"/"+P.locus[m2]->allele2) << " "
       << setw(6) << (P.locus[m2]->allele2+"/"+P.locus[m2]->allele2) << " "
       << setw(6) << "0/0" << " "
       << setw(6) << "*/*" << "\n";

  TWOL << setw(par::pp_maxsnp) << par::twolocus_snp1 << " "
       << setw(4) << (P.locus[m1]->allele1+"/"+P.locus[m1]->allele1) << " ";

  TWOL << setw(6) << count[0] << " "
       << setw(6) << count[1] << " "
       << setw(6) << count[2] << " "
       << setw(6) << count[3] << " "
       << setw(6) << margA[0] << "\n";

  TWOL << setw(par::pp_maxsnp) << " " << " "
       << setw(4) << (P.locus[m1]->allele1+"/"+P.locus[m1]->allele2) << " ";

  TWOL << setw(6) << count[4] << " "
       << setw(6) << count[5] << " "
       << setw(6) << count[6] << " "
       << setw(6) << count[7] << " "
       << setw(6) << margA[1] << "\n";

  TWOL << setw(par::pp_maxsnp) << " " << " "
       << setw(4) << (P.locus[m1]->allele2+"/"+P.locus[m1]->allele2) << " ";

  TWOL << setw(6) << count[8] << " "
       << setw(6) << count[9] << " "
       << setw(6) << count[10] << " "
       << setw(6) << count[11] << " "
       << setw(6) << margA[2] << "\n";

  TWOL << setw(par::pp_maxsnp) << " " << " "
       << setw(4) << "0/0" << " ";

  TWOL << setw(6) << count[12] << " "
       << setw(6) << count[13] << " "
       << setw(6) << count[14] << " "
       << setw(6) << count[15] << " "
       << setw(6) << margA[3] << "\n";

  TWOL << setw(par::pp_maxsnp) << " " << " "
       << setw(4) << "*/*" << " ";

  TWOL << setw(6) << margB[0] << " "
       << setw(6) << margB[1] << " "
       << setw(6) << margB[2] << " "
       << setw(6) << margB[3] << " "
       << setw(6) << total << "\n";

  TWOL << "\n";

}

void Plink::display_twolocus()
{
  
  if (!par::SNP_major)
    Ind2SNP();

  // Find the two SNPs
  
  // Write to file a contingency table, possibly stratified 
  // by phenotype

  string f = par::output_file_name + ".twolocus";
  
  printLOG("Writing two-locus table for " 
	   + par::twolocus_snp1 + " by " 
	   + par::twolocus_snp2 + " to [ " + f + " ] \n");

  ofstream TWOL(f.c_str(), ios::out);
  TWOL.clear();
  TWOL.precision(3);

  int m1 = getMarkerNumber(*this,par::twolocus_snp1);
  int m2 = getMarkerNumber(*this,par::twolocus_snp2);
  
  if (m1<0) error("Marker "+par::twolocus_snp1+" not found\n");
  if (m2<0) error("Marker "+par::twolocus_snp2+" not found\n");
    
  vector<CSNP*>::iterator sa = SNP.begin()+m1;
  vector<CSNP*>::iterator sb = SNP.begin()+m2;
  vector<Locus*>::iterator loc1 = locus.begin()+m1;
  vector<Locus*>::iterator loc2 = locus.begin()+m2;
  vector<bool>::iterator ia1 = (*sa)->one.begin();
  vector<bool>::iterator ia2 = (*sa)->two.begin();
  vector<bool>::iterator ib1 = (*sb)->one.begin();
  vector<bool>::iterator ib2 = (*sb)->two.begin();

  vector<Individual*>::iterator person = sample.begin();
   
  int gtype = 0;

  vector_t c_all(16);
  vector_t c_case(16);
  vector_t c_control(16);
  
  while ( person != sample.end() )
    {
      
      if ( (!*ia1) && (!*ia2) ) 	
	{
	  if      ( (!*ib1) && (!*ib2) ) gtype=0;
	  else if ( (!*ib1) &&   *ib2  ) gtype=1;
	  else if (   *ib1  &&   *ib2  ) gtype=2;
	  else if (   *ib1  && (!*ib2) ) gtype=3;
	}
      else if ( (!*ia1) &&   *ia2  ) 
	{
	  if      ( (!*ib1) && (!*ib2) ) gtype=4;
	  else if ( (!*ib1) &&   *ib2  ) gtype=5;
	  else if (   *ib1  &&   *ib2  ) gtype=6;
	  else if (   *ib1  && (!*ib2) ) gtype=7;
	}
      else if (   *ia1  &&   *ia2  ) 
	{
	  if      ( (!*ib1) && (!*ib2) ) gtype=8;
	  else if ( (!*ib1) &&   *ib2  ) gtype=9;
	  else if (   *ib1  &&   *ib2  ) gtype=10;
	  else if (   *ib1  && (!*ib2) ) gtype=11;
	}
      else if (   *ia1  && (!*ia2) ) 
	{
	  if      ( (!*ib1) && (!*ib2) ) gtype=12;
	  else if ( (!*ib1) &&   *ib2  ) gtype=13;
	  else if (   *ib1  &&   *ib2  ) gtype=14;
	  else if (   *ib1  && (!*ib2) ) gtype=15;
	}

      c_all[gtype]++;
      if ( ! (*person)->missing )
	{
	  if ( (*person)->aff )
	    c_case[gtype]++;
	  else
	    c_control[gtype]++;
	}
      
      ia1++;
      ia2++;
      ib1++;
      ib2++;
      person++;
      
    }


  TWOL << "\nAll individuals\n===============\n";
  make2LTable(TWOL,*this,m1,m2,c_all,false);

  TWOL.setf(ios::fixed);
  make2LTable(TWOL,*this,m1,m2,c_all,true);
  TWOL.unsetf(ios::fixed);

  TWOL << "\nCases\n=====\n";
  make2LTable(TWOL,*this,m1,m2,c_case,false);
  TWOL.setf(ios::fixed);
  make2LTable(TWOL,*this,m1,m2,c_case,true);
  TWOL.unsetf(ios::fixed);
  TWOL << "\nControls\n========\n";
  make2LTable(TWOL,*this,m1,m2,c_control,false);
  TWOL.setf(ios::fixed);
  make2LTable(TWOL,*this,m1,m2,c_control,true);
  TWOL.unsetf(ios::fixed);

  TWOL << "\n";

  TWOL.close();

}


void Plink::extractExcludeSet(bool exclude)
{
  
  // Make map of locus name with 'l' number
  map<string,int> mlocus;
  vector<bool> del(nl_all);
  for (int l=0;l<nl_all;l++)
    {
      mlocus.insert(make_pair(locus[l]->name,l));
      if (exclude) del[l] = false;  // start off all included
      else del[l] = true;           // start off all excluded
    }
    
  map<string,int>::iterator ilocus;
  
  
  //////////////////////////////////
  // Either extract a certain "GENE"

  if (par::dump_gene)
    {

      // Temporarily read sets -- these will get updated later
      readSet();
      
      bool found_gene = false;
      int c=0;

      for (int s=0; s<setname.size(); s++)
	{

	  if ( setname[s] == par::dump_genename )
	    {
	      // Retain just these SNPs
	      for (int l=0; l<snpset[s].size(); l++)
		del[snpset[s][l]] = false;
	      found_gene = true;
	      c = snpset[s].size();
	      break;
	    }
	  if (found_gene) break;
	}
      
      if (found_gene)
	printLOG(int2str(c) + " SNPs selected from [ " 
		 + par::dump_genename + " ]\n");
      else
	printLOG("WARNING: ** Could not find [ " 
		 + par::dump_genename + " ] in [ " 
		 + par::setfile + " ]\n");

    }

  /////////////////////////////////////
  // OR get a list of SNPs form the command line

  else if ( par::snp_include_from_cl ) 
    {
      NList nl(nl_all);
      vector<int> lst = nl.deparseStringList( par::snp_include_range, & mlocus );
      
      for (int l=0; l<lst.size(); l++)
	del[ lst[l] ] = false;

      printLOG("Extracting " + int2str(lst.size() ) +" SNPs\n");
      
    }


  /////////////////////////////////////
  // OR get a list of ranges from a file
  
  else if ( par::snp_range_list )
    {
      

      string filename = exclude ? par::exclude_file : par::extract_file;
      set<int> slist;

      // Read list of ranges     
      map<string, set<Range> > ranges;
      map<int,set<Range*> > snp2range;

      makeScaffold( *this );

      mapRangesToSNPs( filename , 
		       ranges,		       
		       snp2range );

      map<int,set<Range*> >::iterator i1 = snp2range.begin();
      while ( i1 != snp2range.end() )
	{
	  del[ i1->first ] = exclude ? true : false;
	  slist.insert( i1->first );
	  ++i1;
	}

      if ( exclude )
	printLOG("Excluding " + int2str( slist.size() ) + " SNPs\n");
      else
	printLOG("Extracting " + int2str( slist.size() ) + " SNPs\n");
    } 
 
  /////////////////////////////////////
  // OR get a list of SNPs from a file

  else 
    {
      string filename = par::extract_file;
      
      if (exclude)
	filename = par::exclude_file;
      
      checkFileExists(filename);  
      
      ifstream INFILE(filename.c_str(),ios::in);
      INFILE.clear();
      
      printLOG("Reading list of SNPs ");
      
      if (exclude) 
	printLOG("to exclude [ "
		 + par::exclude_file + " ] ... ");
      else 
	printLOG("to extract [ "
		 + par::extract_file + " ] ... ");
      
      int c=0;
      while (!INFILE.eof())
	{
	  string m;
	  INFILE >> m;
	  if (m=="") continue;
	  
	  ilocus = mlocus.find(m);
	  if (ilocus != mlocus.end())
	    {
	      if (exclude) del[ilocus->second] = true; 
	      else del[ilocus->second] = false; 
	      c++;
	    }
	}
      INFILE.close();
      printLOG(int2str(c)+" read\n");
    }
  
  

  ////////////////////////////////////////
  // Remove selected loci from locus list, 

  deleteSNPs(del);
  


}

void Plink::removeIndividuals(bool keep)
{
  

  //////////////////////////////////
  // Make map of individuals FID/IID

  map<string,int> mperson;
  map<string,int>::iterator iperson;
  vector<bool> del(n);
  for (int i=0;i<n;i++)
    {
      mperson.insert(make_pair(sample[i]->fid+"_"+sample[i]->iid,i));
      if (!keep) del[i] = false;    // start off all included
      else del[i] = true;           // start off all excluded
    }
  

  ///////////////////////////
  // Read list of individuals

  string filename = par::remove_indiv_list;
  if (keep)
    filename = par::keep_indiv_list;

  checkFileExists(filename);  
  ifstream INFILE(filename.c_str(),ios::in);
  INFILE.clear();
   
  printLOG("Reading individuals ");
  if (keep) printLOG("to keep [ "+par::keep_indiv_list + " ] ... ");
  else printLOG( "to remove [ "+par::remove_indiv_list + " ] ... ");

  int c=0;
  while (!INFILE.eof())
    {
      vector<string> s = tokenizeLine( INFILE );

      if ( s.size() == 1 ) 
	error("Problem with line:\n"+s[0]);
      else if (s.size() == 0 )
	continue;
      
      string fid = s[0];
      string iid = s[1];

      iperson = mperson.find(fid+"_"+iid);
      if (iperson != mperson.end())
	{
	  if (!keep) del[iperson->second] = true; 
	  else del[iperson->second] = false; 
	  c++;
	}
    }
  INFILE.close();
  printLOG(int2str(c)+" read\n");
  

  ////////////////////////////////////
  // Remove individuals as appropriate

  int n_removed = deleteIndividuals(del);

  printLOG(int2str(n_removed)+" individuals removed with ");
  if (keep) printLOG("--keep option\n");
  else printLOG("--remove option\n");
  
}



void Plink::keep2SetsForGenome()
{
  
  //////////////////////////////////
  // Make map of individuals FID/IID
  
  map<string,int> mperson;
  map<string,int>::iterator iperson;
  vector<bool> del(n);
  for (int i=0;i<n;i++)
    {
      mperson.insert(make_pair(sample[i]->fid+"_"+sample[i]->iid,i));
      del[i] = true;           // start off all excluded
    }
  

  ///////////////////////////
  // Read list of individuals

  checkFileExists(par::genome_setlist1);
  checkFileExists(par::genome_setlist2);

  ifstream INFILE1(par::genome_setlist1.c_str(),ios::in);
  INFILE1.clear();
   
  ifstream INFILE2(par::genome_setlist2.c_str(),ios::in);
  INFILE2.clear();
  
  int c=0;
  while (!INFILE1.eof())
    {
      string fid,iid;
      INFILE1 >> fid >> iid;
      if (fid=="" || iid=="") continue;
      
      iperson = mperson.find(fid+"_"+iid);
      if (iperson != mperson.end())
	{
	  del[iperson->second] = false; 
	  gset1.insert(sample[iperson->second]);
	  c++;
	}
    }
  INFILE1.close();
  printLOG(int2str(c)+" read from [ "+par::genome_setlist1 +" ]\n");
  
  c=0;
  while (!INFILE2.eof())
    {
      string fid,iid;
      INFILE2 >> fid >> iid;
      if (fid=="" || iid=="") continue;
      
      iperson = mperson.find(fid+"_"+iid);
      if (iperson != mperson.end())
	{
	  del[iperson->second] = false; 
	  gset2.insert(sample[iperson->second]);
	  c++;
	}
    }
  INFILE2.close();
  printLOG(int2str(c)+" read from [ "+par::genome_setlist2 +" ]\n");
    
  
  ////////////////////////////////////
  // Remove individuals as appropriate

  int n_removed = deleteIndividuals(del);

}



void Plink::zeroOnCluster()
{
  
  // Make map of locus name with 'l' number
  map<string,int> mlocus;
  vector<bool> del(nl_all);
  for (int l=0;l<nl_all;l++)
    mlocus.insert(make_pair(locus[l]->name,l));
  
  map<string,int>::iterator ilocus;
  
  
  ///////////////////////////////////////////
  // Get a list of SNPs/clusters from a file

  
  string filename = par::zero_cluster_filename;
  checkFileExists(filename);  
  
  ifstream INFILE(filename.c_str(),ios::in);
  INFILE.clear();

  printLOG("Reading list of SNP/clusters to zero out [ ");
  printLOG(par::zero_cluster_filename + " ]\n");
  
  
  int c=0;
  while (!INFILE.eof())
    {
      string m;
      string k;
      INFILE >> m >> k;
      
      if (m=="") continue;
      
      ilocus = mlocus.find(m);
      if (ilocus != mlocus.end())
	{
	  int l = ilocus->second;
	  map<string,int>::iterator ki = kmap.find(k);

	  if ( ki == kmap.end() )
	    {
	      continue; // cluster does not exist
	      //printLOG("Cluster [ "+k+" ] not found \n");
	    }
	  else
	    {
	      int k2 = ki->second;
	      // Zero out the following
	      for (int i=0; i<n; i++)
		{
		  if ( sample[i]->sol == k2 )
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
		      c++;
		    }
		}
	    }
	}
    }
  INFILE.close();
  printLOG(int2str(c)+" genotypes zeroed\n");
}



void Plink::setObligMissing()
{

  // Temporaily use the cluster (sol) variables to store information
  // about oblig missing (i.e. call this *before* we read any 
  // subsequent cluster variables
  

  printLOG("Reading list of SNP/clusters that are obligatory missing [ ");
  printLOG(par::oblig_missing_filename + " ]\n");

  printLOG("Reading the clusters that define obligatory missingness [ ");
  printLOG(par::oblig_clusters_filename + " ]\n");

  string stored_name = par::include_cluster_filename;
  par::include_cluster_filename = par::oblig_clusters_filename;
  readClusterFile();
  par::include_cluster_filename = stored_name;

  // Make map of locus name with 'l' number

  map<string,int> mlocus;
  vector<bool> del(nl_all);
  for (int l=0;l<nl_all;l++)
    mlocus.insert(make_pair(locus[l]->name,l));
  
  map<string,int>::iterator ilocus;
    

  ///////////////////////////////////////////
  // Get a list of SNPs/clusters from a file
  
  string filename = par::oblig_missing_filename;
  checkFileExists(filename);  
  
  ifstream INFILE(filename.c_str(),ios::in);
  INFILE.clear();

  
  int c=0;
  while (!INFILE.eof())
    {
      string m;
      string k;
      INFILE >> m >> k;
      
      if (m=="") continue;
      
      ilocus = mlocus.find(m);
      if (ilocus != mlocus.end())
	{
	  int l = ilocus->second;
	  map<string,int>::iterator ki = kmap.find(k);

	  if ( ki == kmap.end() )
	    continue; // cluster does not exist	      
	  else
	    {
	      int2 p;
	      p.p1 = l;
	      p.p2 = ki->second;
	      oblig_missing.insert(p);
	      ++c;
	    }
	}
    }
  INFILE.close();
  printLOG(int2str(c)+" SNP/cluster combinations set as obligatory missing\n");
}



void Plink::display_recoded_PEDFILE()
{
  
  string f = par::output_file_name + ".ped";
  
  printLOG("Writing recoded ped file to [ " + f + " ] \n");
  ofstream PED(f.c_str(), ios::out);
  PED.clear();

  string missingCode = par::out_missing_phenotype;
  if ( par::recode_HV )
    missingCode = "0";

  for (int i=0;i<n;i++)
    {
      Individual * person = sample[i];
      PED << person->fid << par::recode_delimit
	  << person->iid << par::recode_delimit
	  << person->pat << par::recode_delimit
	  << person->mat << par::recode_delimit
	  << person->sexcode << par::recode_delimit;

      if ( person->missing )
	{
	  PED << missingCode;
	}
      else
	{
	  if (par::bt)
	    PED << (int)person->phenotype;
	  else
	    PED << person->phenotype;
	}

      for (int l=0;l<nl_all;l++)
	PED << genotypeToFile( *this , i, l );

      PED << "\n";
    }
  
  PED.close();
  
  
  // And we also need a new map file
  
  if (par::recode_HV)
    {
      f = par::output_file_name + ".info";
      printLOG("Writing HaploView-format map file [ " + f + " ] \n");
    }
  else
    {
      f = par::output_file_name + ".map";
      printLOG( "Writing new map file to [ " + f + " ] \n");
    }

  ofstream MAP(f.c_str(), ios::out);
  MAP.clear();
  
  bool HV_okay = true;

  int HV_chr = locus[0]->chr;

  if ( par::recode_HV || par::recode_whap )
    {
      if ( nl_all > 5000 ) 
	printLOG(" *** WARNING :  you are exporting a large number of SNPs (>5000) for\n"
		 "                input into Haploview/WHAP -- make sure you have adequate\n"
		 "                system resources to handle this file\n\n");
    }
  
  if ( par::recode_HV)
    for (int l=0;l<nl_all;l++)
      {
	if ( locus[l]->chr != HV_chr ) 
	  HV_okay = false;
	MAP << locus[l]->name << "\t"
	    << locus[l]->bp << "\n";
      }  
  else if ( par::recode_whap )
    {
      for (int l=0;l<nl_all;l++)
	{
	  MAP << locus[l]->chr << "\t"
	      << locus[l]->name << "\t"
	      << locus[l]->bp << "\n";
	}
    }
  else
    for (int l=0;l<nl_all;l++)
      {
	MAP << locus[l]->chr << "\t"
	    << locus[l]->name << "\t"
	    << locus[l]->pos << "\t"
	    << locus[l]->bp << "\n";
      }
  

  MAP.close();

  
  if ( par::recode_HV )
    {
      if ( ! HV_okay )
	printLOG(" *** WARNING :  you've created a Haploview file containing\n"
		 " ***            SNPs from more than 1 chromosome -- these will\n"
		 " ***            not be read properly in Haploview\n\n");
    }


  if ( par::recode_whap ) 
    {
      // WHAP also needs a DAT file
      
      string f = par::output_file_name + ".dat";
      printLOG("Writing whap-format DAT file [ " + f + " ] \n");
      ofstream DAT(f.c_str(), ios::out);
      DAT.clear();

      if ( par::qt ) 
	DAT << "T trait\n";
      else if ( par::coding01 )
	DAT << "B trait\n";
      else
	DAT << "A trait\n";
      
      for (int l=0;l<nl_all;l++)
	DAT << "M " << locus[l]->name << "\n";
      
      DAT.close();
            
    }

  
}


void Plink::display_recoded_MUTLIST()
{
  
  string f = par::output_file_name + ".rlist";
  
  printLOG("Writing rare-list file to [ " + f + " ] \n");

  ofstream PED(f.c_str(), ios::out);
  
  for (int l=0;l<nl_all;l++)
    {
 
      // Hets, homs, missing
      
      set<int> het;
      set<int> hom;
      set<int> missing;
      
      for (int i=0;i<n;i++)
	{
	  Individual * person = sample[i];
	  
	  // A non-reference genotype?
	  
	  bool a1 = par::SNP_major ? SNP[l]->one[i] : person->one[l];
	  bool a2 = par::SNP_major ? SNP[l]->two[i] : person->two[l];
	  
	  if ( (!a1) && (!a2) )
	    hom.insert(i);
	  else if ( (!a1) && a2 )
	    het.insert(i);
	  else if ( a1 && ! a2 ) 
	    missing.insert(i);
	}
      

      // Hets
      if ( het.size() > 0 )
	{
	  PED << locus[l]->name << par::recode_delimit 
	      << "HET" << par::recode_delimit 
	      << locus[l]->allele1 << par::recode_indelimit 
	      << locus[l]->allele2;
	  set<int>::iterator j = het.begin();
	  while ( j != het.end() )
	    {
	      PED << par::recode_delimit << sample[*j]->fid 
		  << par::recode_delimit << sample[*j]->iid;
	      ++j;
	    }
	  PED << "\n";
	}

      // Homozygotes
      if ( hom.size() > 0 )
	{
	  PED << locus[l]->name << par::recode_delimit 
	      << "HOM" << par::recode_delimit 
	      << locus[l]->allele1 << par::recode_indelimit 
	      << locus[l]->allele1;
	  set<int>::iterator j = hom.begin();
	  while ( j != hom.end() )
	    {
	      PED << par::recode_delimit << sample[*j]->fid 
		  << par::recode_delimit << sample[*j]->iid;
	      ++j;
	    }
	  PED << "\n";
	}

      // Missing genotypes
      if ( missing.size() > 0 )
	{
	  PED << locus[l]->name << par::recode_delimit 
	      << "NIL" << par::recode_delimit 
	      << par::missing_genotype << par::recode_indelimit 
	      << par::missing_genotype;
	  set<int>::iterator j = missing.begin();
	  while ( j != missing.end() )
	    {
	      PED << par::recode_delimit << sample[*j]->fid 
		  << par::recode_delimit << sample[*j]->iid;
	      ++j;
	    }
	  PED << "\n";
	}


      }
  
  PED.close();
    
  // And a corresponding MAP file

  f = par::output_file_name + ".map";
  printLOG( "Writing new map file to [ " + f + " ] \n");  
  ofstream MAP(f.c_str(), ios::out);
  for (int l=0;l<nl_all;l++)
    {
      MAP << locus[l]->chr << "\t"
	  << locus[l]->name << "\t"
	  << locus[l]->pos << "\t"
	  << locus[l]->bp << "\n";
    }

  MAP.close();


  printLOG( "Writing pedigree information to [ " 
	    + par::output_file_name + ".fam ] \n");
  ofstream FAM((par::output_file_name+".fam").c_str(), ios::out);

  for (int i=0;i<n;i++)
    {
      Individual * person = sample[i];
      FAM << person->fid << " "
	  << person->iid << " "
	  << person->pat << " "
	  << person->mat << " "
	  << person->sexcode << " ";
      
      if ( person->missing )
	FAM << par::out_missing_phenotype << "\n";
      else
	{
	  if (par::bt)
	    FAM << (int)person->phenotype << "\n";
	  else
	    FAM << person->phenotype << "\n";
	}
    }
  FAM.clear();
  FAM.close();

}

void Plink::display_recoded_LONG()
{
  
  string f = par::output_file_name + ".lgen";
  
  printLOG("Writing recoded LGEN file to [ " + f + " ] \n");

  ofstream PED(f.c_str(), ios::out);

  string missingCode = par::out_missing_phenotype;

  for (int l=0;l<nl_all;l++)
    for (int i=0;i<n;i++)
      {

	Individual * person = sample[i];	

	if ( par::recode_long_ref )
	  {
	    // Skip reference homozygotes?
	    bool a1 = par::SNP_major ? SNP[l]->one[i] : person->one[l];
	    bool a2 = par::SNP_major ? SNP[l]->two[i] : person->two[l];
	    if ( a1 && a2 ) 
	      continue;
	  }
		
	PED << person->fid << par::recode_delimit
	    << person->iid << par::recode_delimit
	    << locus[l]->name << par::recode_delimit
	    << genotypeToFile( *this , i, l ) << "\n";

      }
  
  PED.close();
    
  // And a corresponding MAP file

  f = par::output_file_name + ".map";
  printLOG( "Writing new map file to [ " + f + " ] \n");  
  ofstream MAP(f.c_str(), ios::out);
  for (int l=0;l<nl_all;l++)
    {
      MAP << locus[l]->chr << "\t"
	  << locus[l]->name << "\t"
	  << locus[l]->pos << "\t"
	  << locus[l]->bp << "\n";
    }

  MAP.close();


  //////////////////////////////////////////
  // And a corresponding reference file?
 
  if ( par::recode_long_ref )
    {
      f = par::output_file_name + ".ref";
      printLOG( "Writing new reference file to [ " + f + " ] \n");  
      ofstream MAP(f.c_str(), ios::out);
      for (int l=0;l<nl_all;l++)
	{

	  if ( locus[l]->allele1 == par::missing_genotype && 
	       locus[l]->allele2 == par::missing_genotype )
	    continue;

	  MAP << locus[l]->name;

	  if ( locus[l]->allele2 != par::missing_genotype )
	    MAP << " " << locus[l]->allele2;

	  if ( locus[l]->allele1 != par::missing_genotype )
	    MAP << " " << locus[l]->allele1;
	  MAP << "\n";

	}      
      MAP.close();
    }
  

  //////////////////////////////////////////
  // And a corresponding FAM file?
 
  printLOG( "Writing pedigree information to [ " 
	    + par::output_file_name + ".fam ] \n");
  ofstream FAM((par::output_file_name+".fam").c_str(), ios::out);

  for (int i=0;i<n;i++)
    {
      Individual * person = sample[i];
      FAM << person->fid << " "
	  << person->iid << " "
	  << person->pat << " "
	  << person->mat << " "
	  << person->sexcode << " ";
      
      if ( person->missing )
	FAM << par::out_missing_phenotype << "\n";
      else
	{
	  if (par::bt)
	    FAM << (int)person->phenotype << "\n";
	  else
	    FAM << person->phenotype << "\n";
	}
    }
  FAM.clear();
  FAM.close();
  
}




void Plink::output_fastphase_format()
{

  string f = par::output_file_name + ".recode.phase.inp";
  
  printLOG("Writing fastphase-format file to [ " + f + " ] \n");
  ofstream OUT(f.c_str(), ios::out);
  OUT.clear();

  // Sample size and number of SNPs
  
  OUT << n << "\n" << nl_all << "\nP ";

  // Positions
  
  bool output_okay = true;
  int output_chr = locus[0]->chr;
  
  for (int l=0;l<nl_all;l++)
    {
      if ( locus[l]->chr != output_chr ) 
	output_okay = false;
      OUT << locus[l]->bp << " ";
    }  
  OUT << "\n";

  if ( ! output_okay )
    printLOG(" *** WARNING :  you've created a fastphase file containing\n"
	     " ***            SNPs from more than 1 chromosome -- these will\n"
	     " ***            not be read properly in fastphase\n\n");
  

  for (int i=0;i<n;i++)
    {

      Individual * person = sample[i];

      OUT << "# ID " << person->iid << "\n";
      
      for (int l=0;l<nl_all;l++)
	{

	  bool a1 = par::SNP_major ? SNP[l]->one[i] : person->one[l];
	  bool a2 = par::SNP_major ? SNP[l]->two[i] : person->two[l];

	  if (  ! a1 )  
	    OUT << "0";
	  else 
	    {
	      if ( a2 ) 
		OUT << "1";
	      else
		OUT << "?";
	    }
	}

      OUT << "\n";
      
      // Second allele

      for (int l=0;l<nl_all;l++)
	{
	  
	  bool a1 = par::SNP_major ? SNP[l]->one[i] : person->one[l];
	  bool a2 = par::SNP_major ? SNP[l]->two[i] : person->two[l];

	  if (  ! a1  )  
	    {
	      if ( a2 ) 
		OUT << "1";
	      else 
		OUT << "0";
	    }
	  else 
	    {
	      if ( a2 ) 
		OUT << "1";
	      else
		OUT << "?";
	    }

	}
    
      // Next individual

      OUT << "\n";

    }
    
  OUT.close();

}


void Plink::output_bimbam_format()
{ 
  // One file of SNP locations (SNP BP) *.pos.txt
  // One file of phenotypes (column 6 from FAM)
  // One file containing individual IDs and genotypes (transposed format)
  //   N_ind N_snp
  //   "IND,"IID
  //   'SNP Name' genotypes... 


  // Position file

  string f = par::output_file_name + ".recode.pos.txt";
  printLOG("Writing BIMBAM position file to [ " + f + " ] \n");
  ofstream POS(f.c_str(), ios::out);
  for (int l=0;l<nl_all;l++)
    POS << locus[l]->name << " " << locus[l]->bp << "\n";
  POS.close();


  // Phenotype information

  f = par::output_file_name + ".recode.pheno.txt";
  printLOG("Writing BIMBAM phenotype file to [ " + f + " ] \n");
  ofstream PHE(f.c_str(), ios::out);

  if (par::bt)
    for (int i=0;i<n;i++)
      PHE << (int)sample[i]->phenotype << "\n";
  else
    for (int i=0;i<n;i++)
      PHE << sample[i]->phenotype << "\n";
  
  PHE.close();
  

  // Genotype file

  f = par::output_file_name + ".recode.geno.txt";
  printLOG("Writing BIMBAM genotype file to [ " + f + " ] \n");
  ofstream GEN(f.c_str(), ios::out);

  par::missing_genotype = "?";
  par::recode_delimit = ",";
  par::recode_indelimit = "";

  GEN << n << "\n"
      << nl_all << "\n"
      << "IND";
  
  for (int i=0;i<n;i++)
    GEN << "," << sample[i]->iid;
  
  GEN << "\n";

  for (int l=0;l<nl_all;l++)
    {
      GEN << locus[l]->name;
      
      for (int i=0;i<n;i++)
	GEN << genotypeToFile( *this , i, l);

      GEN << "\n";
    }
  
  GEN.close();

}

void Plink::output_structure_format()
{

  // SNP names; inter SNP distances (with -1 starting each new
  // chromosome); 


  string f = par::output_file_name + ".recode.strct_in";
  printLOG("Writing STRUCTURE format file to [ " + f + " ] \n");
  ofstream GEN(f.c_str(), ios::out);

  for (int l=0;l<nl_all;l++)
    GEN << locus[l]->name << " ";
  GEN << "\n";

  for (int l=0;l<nl_all;l++)
    GEN << 
      ( l==0 || locus[l-1]->chr != locus[l]->chr ? 
	-1 :
	locus[l]->bp - locus[l-1]->bp ) 
	<< " ";
  GEN << "\n";

  par::missing_genotype = "0";
  par::recode_delimit = " ";
  par::recode_indelimit = " ";
  par::recode_12 = true;

  map<string,int> fmap;
  int cnt = 0;
  for (int i=0;i<n;i++)
    {
      Individual * person = sample[i];
      if ( fmap.find( person->fid ) == fmap.end() )
	fmap.insert(make_pair(person->fid,++cnt));
    }

  for (int i=0;i<n;i++)
    {
      GEN << sample[i]->iid << " "
	  << fmap.find(sample[i]->fid)->second;
      
      for (int l=0;l<nl_all;l++)
	GEN << genotypeToFile( *this , i, l);
      
      GEN << "\n";
    }
  
  GEN.close();


//   GEN << n << "\n"
//       << nl_all << "\n"
//       << "IND";
  
}


void Plink::display_recoded_PEDFILE_transpose()
{

  string f = par::output_file_name + ".tped";
  printLOG("Writing transposed ped file to [ " + f + " ] \n");
  ofstream PED(f.c_str(), ios::out);
  PED.clear();

  //  par::recode_indelimit = "";

  for (int l=0;l<nl_all;l++)
    {
      PED << locus[l]->chr << par::recode_delimit
 	  << locus[l]->name << par::recode_delimit
 	  << locus[l]->pos << par::recode_delimit
	  << locus[l]->bp;
      
      for (int i=0;i<n;i++)
	PED << genotypeToFile( *this , i, l);

      PED << "\n";
    }
  
  PED.close();
  
  
  // And we also need a new family info file

  f = par::output_file_name + ".tfam";

  printLOG( "Writing family information to [ " + f + " ] \n");  
  ofstream FAM(f.c_str(), ios::out);
  FAM.clear();
  
  for (int i=0;i<n;i++)
    {
      Individual * person = sample[i];

      FAM << person->fid << par::recode_delimit
	  << person->iid << par::recode_delimit
	  << person->pat << par::recode_delimit
	  << person->mat << par::recode_delimit
	  << person->sexcode << par::recode_delimit;
      
      if ( person->missing )
 	FAM << par::out_missing_phenotype;
      else
	{
	  if (par::bt)
	    FAM << (int)person->phenotype;
	  else
	    FAM << person->phenotype;
	}
      
      FAM << "\n";
    }
  
  FAM.close();
  
}

void Plink::display_recoded_PEDFILE_AD()
{

  string f = par::output_file_name + ".raw";
  printLOG( "Writing recoded file to [ " + f + " ] \n");
  ofstream PED(f.c_str(), ios::out);
  PED.clear();

  map<string,string> amap;
  if ( par::recode_allele_coding )
    {
      string f = par::recode_allele_coding_file;
      printLOG( "Reading allele coding list from [ " + f + " ] \n");
      ifstream AMAP;
      checkFileExists(f);
      AMAP.open(f.c_str(), ios::in);
      while ( ! AMAP.eof() )
	{
	  string snp, allele;
	  AMAP >> snp >> allele;
	  if (snp=="") break;
	  amap.insert(make_pair(snp,allele));
	}
      printLOG("Read allele codes for " + int2str( amap.size() ) + " SNPs\n");
      AMAP.close();
    }

  // Header row
  PED << "FID" 
      << par::recode_delimit << "IID"
      << par::recode_delimit << "PAT"
      << par::recode_delimit << "MAT"
      << par::recode_delimit << "SEX"
      << par::recode_delimit << "PHENOTYPE";

  for (int l=0;l<nl_all;l++)
    {
      
      string aname = locus[l]->allele1;

      if ( par::recode_allele_coding )
	{
	  map<string,string>::iterator a = amap.find(locus[l]->name);
	  if ( a != amap.end() )
	    aname = a->second;
	}

      PED << par::recode_delimit << locus[l]->name+"_"+aname;
      
      if ( ! par::recode_AD_Aonly )
	PED << par::recode_delimit << locus[l]->name+"_HET";
    }
  
  PED << "\n";
  
  for (int i=0;i<n;i++)
    {
      Individual * person = sample[i];
      PED << person->fid 
	  << par::recode_delimit << person->iid 
	  << par::recode_delimit << person->pat 
	  << par::recode_delimit << person->mat 
	  << par::recode_delimit << person->sexcode;

      if ( person->missing )
	PED << par::recode_delimit << par::out_missing_phenotype;
      else
	{
	  if (par::bt)
	    PED << par::recode_delimit << (int)person->phenotype;
	  else
	    PED << par::recode_delimit << person->phenotype;
	}

      string g0 = par::recode_AD_Aonly ? 
	par::recode_delimit + "2" 
	: par::recode_delimit + "2" + par::recode_delimit + "0";
      
      string g1 = par::recode_AD_Aonly ? 
	par::recode_delimit + "1" 
	: par::recode_delimit + "1" + par::recode_delimit + "1";
      
      string g2 = par::recode_AD_Aonly ? 
	par::recode_delimit + "0" 
	: par::recode_delimit + "0" + par::recode_delimit + "0";
	  
      string gX = par::recode_AD_Aonly ? 
	par::recode_delimit + "NA" 
	: par::recode_delimit + "NA" + par::recode_delimit + "NA";
      
      
      for (int l=0;l<nl_all;l++)
	{
	  
	  bool a1 = par::SNP_major ? 
	    SNP[l]->one[i] : 
	    sample[i]->one[l];
	  
	  bool a2 = par::SNP_major ? 
	    SNP[l]->two[i] : 
	    sample[i]->two[l];
	  
	  if (par::recode_AD_fixed && locus[l]->allele1=="1")
	    {
	      if ( (!a1) && (!a2) )  
		PED << g2;
	      else if ( (!a1) && a2) 
		PED << g1;
	      else if ( a1 && a2 ) 
		PED << g0;
	      else 
		PED << gX;
	    }
	  else if ( par::recode_allele_coding )
	    {
	      Locus * loc = locus[l];
	      map<string,string>::iterator a = amap.find(loc->name);
	      bool flip = false;

	      if ( a != amap.end() )
		{
		  if ( a->second == loc->allele2 )
		    flip = true;
		}
	      
	      if ( flip ) 
		{
		  if ( (!a1) && (!a2) )  
		    PED << g2;
		  else if ( (!a1) && a2) 
		    PED << g1;
		  else if ( a1 && a2 ) 
		    PED << g0;
		  else 
		    PED << gX;	      
		}
	      else
		{
		  if ( (!a1) && (!a2) )  
		    PED << g0;
		  else if ( (!a1) && a2) 
		    PED << g1;
		  else if ( a1 && a2 ) 
		    PED << g2;
		  else 
		    PED << gX;	      
		}

	    }
	  else  
	    {
	      if ( (!a1) && (!a2) )  
		PED << g0;
	      else if ( (!a1) && a2) 
		PED << g1;
	      else if ( a1 && a2 ) 
		PED << g2;
	      else 
		PED << gX;	      
	    }
	}

      
      PED << "\n";
    }
  
  PED.close();

}



void Plink::write_covariates()
{
  
  printLOG( "Writing covariate information to [ " 
	    + par::output_file_name + ".cov ] \n");
  ofstream COV((par::output_file_name+".cov").c_str(), ios::out);
  
  // First, make dummy-variables for any multi-category variables; up
  // to a limit (10 levels)

  vector<int> downcoding_level(par::clist_number);
  vector< map<int,int> > levels(par::clist_number);
  vector< map<int,int> > backcode(par::clist_number);
  
  if ( par::dump_covar_dummy_coding ) 
    {
      for (int c=0; c<par::clist_number; c++)
	{	  
	  
     	  // only look at non-missing individuals (should be all okay
     	  // for covariate)

	  for (int i=0;i<n;i++)
	    {
	      
	      if ( ! sample[i]->clistMissing[c] )
		{
		  int covariate = (int)sample[i]->clist[c];
		  
		  if ( levels[c].find( covariate ) == levels[c].end() )
		    {
		      int t = levels[c].size();
		      levels[c].insert( make_pair(covariate, t));
		      backcode[c].insert( make_pair(t,covariate));
		    }
		}
	    }
	  
	  if ( levels[c].size() > 2 && levels[c].size() < 50 ) 
	    downcoding_level[c] = levels[c].size();
	  else
	    downcoding_level[c] = 0;
	  
	}
      
    }
  

  // Always write a header row
  
  COV << "FID IID ";
  if (par::dump_covar_with_phenotype)
    COV << "PAT MAT SEX PHENOTYPE ";  

  // Covariate names in header row (possibly with dummy-variable
  // downcoding)

  if ( par::dump_covar_dummy_coding )
    {
      for (int c=0; c<par::clist_number; c++)
	{
	  if ( downcoding_level[c] == 0 )
	    COV << clistname[c] << " ";
	  else
	    {
	      for (int d=1; d<downcoding_level[c]; d++) // skip first
		COV << clistname[c] << "_" << backcode[c].find(d)->second << " ";	      
	    }
	}
    }
  else
    for (int c=0; c<par::clist_number; c++)
      COV << clistname[c] << " ";

  COV << "\n";
  
  
  // Output values for each individual

  for (int i=0;i<n;i++)
    {
      Individual * person = sample[i];
      
      COV << person->fid << " "
	  << person->iid << " ";

      if (par::dump_covar_with_phenotype)
	{
	  COV << person->pat << " "
	      << person->mat << " "
	      << person->sex << " ";
	  
	  if ( person->missing )
	    COV << par::out_missing_phenotype << " ";
	  else
	    COV << person->phenotype << " ";
	}
      
      if ( par::dump_covar_dummy_coding )
	{
	  for (int c=0; c<par::clist_number; c++)
	    {
	      if ( downcoding_level[c] == 0 )
		{
		  if ( person->clistMissing[c] )
		    COV << par::out_missing_phenotype << " ";
		  else      
		    COV << person->clist[c] << " ";		  
		}
	      else
		{
		  for (int d=1; d<downcoding_level[c]; d++)
		    {
		      if ( person->clistMissing[c] )
			COV << par::out_missing_phenotype << " ";
		      else
			{
			  if ( levels[c].find( (int)person->clist[c] )->second == d ) 
			    COV << "1 ";
			  else
			    COV << "0 ";
			}
		    }
		}
	    }
	}
      else
	{
	  for (int c=0; c<par::clist_number; c++)
	    COV << person->clist[c] << " ";
	}
      COV << "\n";
    }

  COV.close();
}


void Plink::write_clusters()
{
  
  printLOG( "Writing cluster information to [ " 
	    + par::output_file_name + ".clst ] \n");
  ofstream COV((par::output_file_name+".clst").c_str(), ios::out);
  


  // Do not write a header row
  
  //  COV << "FID IID CLST";
  // COV << "\n";
  
  
  // Output values for each individual

  for (int i=0;i<n;i++)
    {
      Individual * person = sample[i];
      
      COV << person->fid << " "
	  << person->iid << " ";
      if ( person->sol == -1 )
	COV << "NA\n";
      else
	COV << kname[ person->sol ] << "\n";
    }

  COV.close();
}


void Plink::write_snplist()
{
  printLOG( "Writing list of SNPs to [ " 
	    + par::output_file_name + ".snplist ] \n");
  ofstream SL((par::output_file_name+".snplist").c_str(), ios::out);
  for (int l=0;l<nl_all;l++)
    SL << locus[l]->name << "\n";
  SL.close();
}


void Plink::write_BITFILE()
{

  printLOG( "Writing pedigree information to [ " 
	    + par::output_file_name + ".fam ] \n");
  ofstream BIT((par::output_file_name+".fam").c_str(), ios::out);
  // For each individual

  for (int i=0;i<n;i++)
    {
      Individual * person = sample[i];
      BIT << person->fid << " "
	  << person->iid << " "
	  << person->pat << " "
	  << person->mat << " "
	  << person->sexcode << " ";
      
      if ( person->missing )
	BIT << par::out_missing_phenotype << "\n";
      else
	{
	  if (par::bt)
	    BIT << (int)person->phenotype << "\n";
	  else
	    BIT << person->phenotype << "\n";
	}
    }
  BIT.clear();
  BIT.close();
  
  printLOG( "Writing map (extended format) information to [ " 
	    + par::output_file_name + ".bim ] \n");
  BIT.open((par::output_file_name+".bim").c_str(), ios::out);
  
  for (int l=0;l<nl_all;l++)
    {

      if (locus[l]->allele1=="")
 	{
 	  if (locus[l]->allele2!="0") locus[l]->allele1="0";
 	  else locus[l]->allele1="X";
 	}
      
      if (locus[l]->allele2=="")
	locus[l]->allele2="0";
      
      BIT << locus[l]->chr << "\t"
	  << locus[l]->name << "\t"
	  << locus[l]->pos << "\t"
	  << locus[l]->bp << "\t" 
	  << locus[l]->allele1 << "\t"
	  << locus[l]->allele2 << "\n";
    }
  
  BIT.clear();
  BIT.close();


  //////////////////////////////////////
  // Save genotype data in BITFILE format

  printLOG("Writing genotype bitfile to [ " 
	   + par::output_file_name + ".bed ] \n");
    
  BIT.open((par::output_file_name+".bed").c_str(), ios::out | ios::binary);
  
  if (par::out_SNP_major)
    printLOG("Using (default) SNP-major mode\n");
  else
    printLOG("Using individual-major mode\n");

  bitset<8> b;
  char ch[1];


  // Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file

  b.reset();
  b.set(2);  b.set(3);  b.set(5);  b.set(6);
  ch[0] = (char)b.to_ulong();
  BIT.write(ch,1);

  b.reset();
  b.set(0);  b.set(1);  b.set(3);  b.set(4);
  ch[0] = (char)b.to_ulong();
  BIT.write(ch,1);


  // BIT represents status of SNP-major (true) or Ind-major (false)

  b.reset();
  if (par::out_SNP_major)
    b.set(0);
  ch[0] = (char)b.to_ulong();
  BIT.write(ch,1);
  
  
  //////////////////////////////////////////
  // Now consider genotypes: SNP-major mode
  
  if (par::out_SNP_major)
    {  
      
      // Create the SNP-major dataset, if need be
      if (!par::SNP_major) Ind2SNP();

      vector<CSNP*>::iterator s = SNP.begin();
         
      // Outer loop over SNPs
      while ( s != SNP.end() )
	{
	  
	  vector<bool>::iterator i1 = (*s)->one.begin();
	  vector<bool>::iterator i2 = (*s)->two.begin();
	  //	  vector<Individual*>::iterator person = sample.begin();

	  // Inner loop over individuals
	  while ( i1 != (*s)->one.end() )
	    {
	      bitset<8> b;
	      b.reset();
	      int c=0;      
	      
	      while (c<8 && i1 != (*s)->one.end() )
		{
		  if ( *(i1) ) b.set(c);
		  i1++;
		  c++;

		  if ( *(i2) ) b.set(c);		  
		  i2++;
		  c++;

		  //person++;
		}
	      
	      char ch[1];
	      ch[0] = (char)b.to_ulong();
	      BIT.write(ch,1);
	      
	    }

	  // next SNP
	  s++;
	}
    }
  

  ////////////////////////////////////////////////////////////
  // Alternatively, consider genotypes: Individual-major mode

  else
    {

      // Create the individual-major dataset, if need be
      if (par::SNP_major) SNP2Ind();
      
      vector<Individual*>::iterator person = sample.begin();
            
      // Outer loop over individuals
      while ( person != sample.end() )
	{
	  
	  vector<bool>::iterator i1 = (*person)->one.begin();
	  vector<bool>::iterator i2 = (*person)->two.begin();
	  
	  // Inner loop over SNPs
	  while ( i1 != (*person)->one.end() )
	    {
	      
	      bitset<8> b;
	      b.reset();
	      int c=0;      
	      
	      while (c<8 && i1 != (*person)->one.end() )
		{
		  if ( *(i1) ) b.set(c);
		  i1++;
		  c++;
		  if ( *(i2) ) b.set(c);		  
		  i2++;
		  c++;
		}
	      
	      char ch[1];
	      ch[0] = (char)b.to_ulong();
	      BIT.write(ch,1);
	      
	    }
	  
	  // next person
	  person++;
	}
      

    }


  BIT.close();
  
}


void Plink::outputSetFile()
{
  
  
  if ( par::make_set_complement )
    {
      if ( !par::make_set_ignore_group )
	printLOG("Making sets of SNPs not in each group\n");
      else
	printLOG("Making a set of all SNPs not in region(s)\n");
    }
  else if ( par::make_set_collapse )
    {
      if ( !par::make_set_ignore_group )
	printLOG("Collapsing all regions into groups\n");
      else
	printLOG("Collapsing all regions into a single set\n");

    }
  

  // We need a scaffold in place now, for range lookup

  makeScaffold(*this);

  map<string, set<Range> > ranges = readRange( par::make_set_file );
  map<string, set<Range> >::iterator r = ranges.begin();
  
  int eset = 0 , fset = 0;

  // Clear set storage
  
  for (int i=0; i<snpset.size();i++)
    snpset[i].clear();
  snpset.clear();
  setname.clear();
  vector<set<int> > inSet;
  
  if ( par::make_set_complement || par::make_set_collapse )
    {

      // Collapse to a single group...

      if ( par::make_set_ignore_group )
	{
	  setname.resize(1);
	  inSet.resize(1);
	  setname[0] = par::make_set_collapse_label;
	}
      else
	{

	  // ...or collapse to "group" label level?

	  setname.resize( Range::groupNames.size() );
	  inSet.resize( Range::groupNames.size() );
	  map<string,int>::iterator i1 = Range::groupNames.begin();
	  while ( i1 != Range::groupNames.end() )
	    {
	      setname[i1->second] = par::make_set_collapse ? 
		i1->first : "C_" + i1->first;
	      ++i1;
	    }
	}
    }
  


  // Loop over all ranges

  int cnt = -1;
  while ( r != ranges.end() )
    {
      
      
      set<Range> * theseRanges = &( r->second );

      // Assume a single, unqiue named range (and so set-size=1)
      
      set<Range>::iterator thisRangeSet = theseRanges->begin();

      const Range * thisRange = &(*thisRangeSet);

      // Get group name/label

      if ( ! ( par::make_set_complement || par::make_set_collapse ) )
	{
	  set<int> t;
	  inSet.push_back(t);
	  ++cnt;
	  setname.push_back( thisRange->name );
	}
      else 
	{
	  // Collapsing, but by gorup
	  if ( ! par::make_set_ignore_group )
	    {
	      cnt = thisRange->group;
	    }
	}
      

      // Consider all inner ranges in this set
      
      while ( thisRangeSet != theseRanges->end() )
	{
	  
	  const Range * thisRange = &(*thisRangeSet);
	  
	  // Assign SNPs
	  int2 srange = mapSNPs2Range(*this, thisRange);
	  
	  
	  // Write SNPs
	  if ( srange.p1 != -1 ) 
	    {
	      
	      if ( par::make_set_complement || par::make_set_collapse ) 
		{
		  if ( par::make_set_ignore_group )
		    for ( int l = srange.p1; l<= srange.p2; l++)
		      inSet[0].insert(l);
		  else
		    for ( int l = srange.p1; l<= srange.p2; l++)
		      inSet[cnt].insert(l);
		}
	      else
		{	      
		  for ( int l = srange.p1; l<= srange.p2; l++)
		    inSet[cnt].insert(l);
		}	  
	      ++fset;
	    }
	  else
	    ++eset; 
	  
	  // Next range in this set
	  ++thisRangeSet;
	}
      // Next set of ranges
      ++r;
    }      
  


  printLOG("Found " + int2str( fset ) 
	   + " ranges with at least 1 SNP; " 
	   + int2str(eset) + " empty ranges\n");


  
  // Copy SNPs into sets

  snpset.resize( inSet.size() );
  for (int j=0; j<inSet.size(); j++)
    {
      snpset[j].clear();

      if ( par::make_set_complement )
	{
	  set<int> & thisSet = inSet[j];
	  for (int l=0; l<nl_all; l++)
	    {
	      set<int>::iterator i = thisSet.find(l);
	      if ( i == thisSet.end() )
		{
		  snpset[j].push_back( l );
		}
	      
	    }
	}
      else
	{
	  set<int>::iterator i = inSet[j].begin();
	  while ( i != inSet[j].end() )
	    {
	      snpset[j].push_back( *i );
	      ++i;
	    }
	}
    }
  
}


void Plink::setTable()
{


  // Write SNPs out, scoring with sets (0/1)
  
  printLOG("Writing set in table form to [ " 
	   + par::output_file_name 
	   + ".set.table ]\n");

  ofstream O( ( par::output_file_name+".set.table" ).c_str() , ios::out );

  O << "SNP\tCHR\tBP";

  for (int i=0; i < setname.size();i++)
    O << "\t" << setname[i]; 
  O << "\n";

  pS->initialiseSetMapping();
  
  for (int l=0; l<nl_all; l++)
    {

      O << locus[l]->name << "\t"
	<< locus[l]->chr << "\t"
	<< locus[l]->bp << "\t";
      
      for (int i=0;i<setname.size();i++)
	{

	  map<int, set<int> >::iterator si = pS->setMapping.find(l);

	  set<int>::iterator si2 = si->second.find(i);

	  if ( si2 == si->second.end() )
	    O << "\t" << "0";
	  else
	    O << "\t" << "1";	  
	}

      O << "\n";
    }
}



void Plink::writeSetFile()
{


  // Write SNPs out, scoring with sets (0/1)
  
  printLOG("Writing set file to [ " + par::output_file_name + ".set ]\n");

  ofstream O( ( par::output_file_name+".set" ).c_str() , ios::out );

  for (int j=0; j<snpset.size(); j++)
    {
      O << setname[j] << "\n";
      for (int l=0; l<snpset[j].size(); l++)
	O << locus[ snpset[j][l] ]->name << "\n";
      O << "END\n\n";
    }

  O.close();

}


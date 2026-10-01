

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
#include <fstream>
#include <sstream>
#include "plink.h"
#include "options.h"
#include "helper.h"

class Pair {
public:
  Individual * p1;
  Individual * p2;

  bool operator< (const Pair & b) const
  {
    return ( ( p1 < b.p1 ) || ( p1 == b.p1 && p2 < b.p2 ) );
  }
  
};

int Plink::readInformative()
{

  skip_pair.clear();


  ////////////////////////////////////
  // Request to read from .genome file
  
  int c=0;
  
  checkFileExists(par::ibd_file);
  printLOG("Reading genome-wide IBD estimates from [ "+par::ibd_file+" ] \n");

  ifstream INC;
  INC.open(par::ibd_file.c_str(), ios::in);

  map<string,Individual*> mperson;
  vector<Individual*>::iterator person = sample.begin();
  while ( person != sample.end() )
    {
      mperson.insert(make_pair( (*person)->fid+"_"+(*person)->iid , *person )); 
      person++;
    }
  
  map<Pair,Z> mpair;

  // Read in .genome file -- from header, get field values for: 
  // FID1 -- IID2 and 
  // z0, z1, z2

  int z0_code = -1;
  int z1_code = -1;
  int z2_code = -1;

  int col_length = 0;
  
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
      if ( tokens[i] == "Z0" )
	z0_code = i;
      if ( tokens[i] == "Z1" )
	z1_code = i;
      if ( tokens[i] == "Z2" )
	z2_code = i;
    }
  
  if ( z0_code == -1 || z1_code == -1 || z2_code == -1 )
    error("Could not find Z0, Z1 or Z2 fields in .genome file");

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

      if (fid1=="") continue;
      

      string z0 = tokens[z0_code];
      string z1 = tokens[z1_code];
      string z2 = tokens[z2_code];

      Z z;
      
      if ( ! ( from_string<double>( z.z0 , z0 , std::dec) && 
	       from_string<double>( z.z1 , z1 , std::dec) && 
	       from_string<double>( z.z2 , z2 , std::dec) ) )	     
	{
	  z.z0 = 1;
	  z.z1 = 0;
	  z.z2 = 0;
	}

      if ( par::debug ) 
	cerr << "Read from file: " 
	     << fid1 << " "
	     << iid1 << ", "
	     << fid2 << " "
	     << iid2 << ", "
	     << z.z0 << " " << z.z1 <<" " << z.z2 << "\n";

	  

      ///////////////////////////////////
      // Range of Genome-Wide IBD okay?
      
      bool val = true;
      double pihat = z.z1/2 + z.z2;
      if ( pihat < par::MIN_PIHAT )
	{
	  if (par::include_all_pairs)
	    {
	      z.z0 = 1-par::include_all_z1;    
	      z.z1 = par::include_all_z1;
	      z.z2 = 0;		      
	    }
	  else
	    val = false;
	}
      
      
      // Above IBD threshold?
      if ( z.z1/2 + z.z2 > par::MAX_PIHAT ) 
	val = false;
      
      // Not a parent-offspring pair?
      if (z.z1 > 0.9 ) val = false;
      
      // Need to nudge?
      if (par::nudge && ( pihat * pihat ) < z.z2 )
	{
	  z.z0 = ( 1 - pihat) * ( 1 - pihat);
	  z.z1 = 2 * pihat * (1-pihat);
	  z.z2 = pihat * pihat;
	}
      
      if (val)  
	{

	  map<string,Individual*>::iterator person1 = mperson.find(fid1+"_"+iid1);
	  map<string,Individual*>::iterator person2 = mperson.find(fid2+"_"+iid2);	  

	  if ( person1 == mperson.end() || person2 == mperson.end() ) 
	    continue;

	  Pair p;
	  p.p1 = person1->second;
	  p.p2 = person2->second;
	  mpair.insert(make_pair(p,z));

	}

    }

  INC.close();



  // Now we've finished reading in all known genome-wide IBD values, 
  // and we've saved those that are in the desired range.
  
  // Consider all pairs
  for (int i1=0; i1<n-1; i1++)
    for (int i2=i1+1; i2<n; i2++)
      {
	
	// Did we see this pair?
	Pair p;
	p.p1 = sample[i1];
	p.p2 = sample[i2];

	map<Pair,Z>::iterator i = mpair.find(p);
	
	// Pair either not found...
	if ( i == mpair.end() )
	  {
	    skip_pair.push_back(true);		    
	  }
	else // ... or was in .genome with valid IBD
	  {
	    skip_pair.push_back(false);
	    c++;
	    saved_IBDg.push_back(i->second);	    	    

	    // And related enough for the --genome-test?

	    if ( par::genome_test ) 
	      {
		Z ibd = i->second;
		if ( ibd.z1/2 + ibd.z2 >= par::genome_test_threshold )
		  {
		    int2 pair;
		    pair.p1 = i1;
		    pair.p2 = i2;
		    related.insert(pair);
		  }
	      }

	  }
      }

  stringstream s2;
  s2 << "\n" << c << " pairs are informative  ( " 
     << par::MIN_PIHAT 
     << " <= pihat <= " 
     << par::MAX_PIHAT        
     << " )\n";
  printLOG(s2.str());
  
  return c;
}


int Plink::calcInformative() 
{
  //////////////////////////////
  // Precalculate to get list of 
  // individuals to skip in subsequent 
  // sections
  
  int c=0;
  int c0=0;

  printLOG("Preprocessing to assess number of informative pairs ... \n");
  for (int i1=0; i1<n-1; i1++)
    for (int i2=i1+1; i2<n; i2++)
      {	    

	c0++;

	Individual * p1 = sample[i1];
	Individual * p2 = sample[i2];

	if (!p1->missing && !p2->missing)
	  {
	    
	    Z IBSg;
	    Z IBDg;
	    
	    // Calculate or fix
	    if (!par::FIXED)
	      {
		IBSg = calcGenomeIBS(p1,p2);
		IBDg = calcGenomeIBD(p1,p2,IBSg);
	      }
	    else
	      {
		IBDg = par::FIX_IBD;
	      }
	    
	    	    
	    bool val = true;

	    //////////////////////////////
	    // Do we meet criteria?
	    
	    // Range of Genome-Wide IBD okay?
	    if (IBDg.z1/2 + IBDg.z2 < par::MIN_PIHAT )
	      {
		if (par::include_all_pairs)
		  {
		    IBDg.z0 = 1-par::include_all_z1;    
		    IBDg.z1 = par::include_all_z1;
		    IBDg.z2 = 0;		      
		  }
		else
		  val = false;
	      }

	    // Above IBD threshold?
	    if ( IBDg.z1/2 + IBDg.z2 > par::MAX_PIHAT ) 
	      val = false;
	    
	    // Not a parent-offspring pair?
	    if (IBDg.z1 > 0.9 ) val = false;
	    
	    // Affected-only pair analysis?
	    
	    if (val)
	      {
		saved_IBDg.push_back(IBDg);
		skip_pair.push_back(false);
		++c;
		cout << c << " pairs extracted; "
		     << c0 << " processed of " << np 
		     << "            \r";
		cout.flush();
	       		
	      }
	    else
	      skip_pair.push_back(true);
	  }
	else
	  skip_pair.push_back(true);
      }		    
  
  
  stringstream s2;
  s2 << "\n" << c << " pairs are informative  ( " 
     << par::MIN_PIHAT 
     << " <= pihat <= " 
     << par::MAX_PIHAT        
     << " )\n";
  printLOG(s2.str());
  
  return c;
}


void Plink::writeInformative()
{
  // Replaced by standard .genome output
}


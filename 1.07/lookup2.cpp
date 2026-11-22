

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
#include <sstream>

#include "plink.h"
#include "helper.h"
#include "options.h"
#include "sockets.h"
#include "nlist.h"

using namespace std;

#define  PORT_NUM                80     
#define  IP_ADDR    "152.19.78.148"  



void convertPosition(string pquery, int & chr, int & bp1, int & bp2, bool useKb, bool useMb)
{
  
  size_t p1 = pquery.find(":");
  size_t p2 = pquery.find("-");
	   
  if ( p1 == string::npos || p2 == string::npos || p1 >= p2 )
    error("Badly formed positional query: " + pquery );
  
  string ccode = pquery.substr(0,p1);
  if ( ccode.size() < 4 )
    error("Badly formed positional query: " + pquery );
  if ( ccode.substr(0,3) != "chr" )
    error("Badly formed positional query: " + pquery );
  
  if ( ! from_string<int>( chr , ccode.substr(3) , std::dec ) )
    error("Badly formed positional query: " + pquery );
  
  string p1code = pquery.substr(p1+1,p2-p1-1);
  double pp1;
  if ( ! from_string<double>( pp1 , p1code , std::dec ) )
    error("Badly formed positional query: " + pquery );
  
  string p2code = pquery.substr(p2+1);
  double pp2;
  if ( ! from_string<double>( pp2 , p2code , std::dec ) )
    error("Badly formed positional query: " + pquery );

  if ( useKb ) 
    {
      bp1 = int(pp1 * 1000);
      bp2 = int(pp2 * 1000);
    }
  else if ( useMb ) 
    {
     bp1 = int(pp1 * 1000000);
     bp2 = int(pp2 * 1000000);
    }
  else
    {
      bp1 = int(pp1);
      bp2 = int(pp2);
    }

  // Possible overflow
  if ( bp1 < 0 ) bp1 = 0;
  if ( bp2 < 0 ) bp2 = 0;
}

void Plink::lookup2()
{

  // In general, for these lookups, we do not want to treat the
  // minus/hyphen character as a range delimiter in the normal
  // sense. It will be specially handled for chr1:1-100 values; also,
  // it might appear in gene names, e.g. HLA-A.  By setting the range
  // character to a space, we essentially ensure that we will never
  // encounter it.  Because we shutdown() after performing this
  // analysis, we do not need to worry about other operations getting
  // messed up.
  
  par::range_delimiter = " ";
    
#ifdef SKIP
   printLOG("Web-lookup not implemented on this system...\n");
   return;
#else

//   printLOG("PLINK-SNP (WGAS SNP annotation courtesy of Patrick Sullivan)\n");
  
   // http://sullivanlab.unc.edu/plink/snp.php?rsid=rs12345,rs67890&gene=DISC1,CACNA1C&exon=1&cnv=0

   // http://sullivanlab.unc.edu/plink/genelist.php?&searchtype=gn&gene=COMT,GTF2H2,FURBERG,AMY1,AMY2B
   // http://sullivanlab.unc.edu/plink/genelist.php?searchtype=gp&chr=6&start=30000000&end=30500000

   //  { rsid = rs12345 OR rsid = rs67890 }
   //  AND { exon = true }
   // AND { cnv = false }
   // AND { gene = DISC1 OR gene = CACNA1C }

   
   string GET_STRING;
   string command;


// Lookup genes

// # 'pos' key -- default is BP
// --lookup-gene [pos,mb] chr18:30.2-30.2 
// --lookup-gene [pos,kb] chr18:30200-302809
// --lookup-gene [pos] chr18:30200000-30280900

// # list

// --lookup-gene [file] mygenes.lst
// --lookup-gene [list] mygenes.lst

// # name

// --lookup-gene CACNA1C
// --lookup-gene [name] CACNA1C

// # query
// --lookup-gene [query] gset=GO_363526,pos=chr6:30200000-30280900,list=mygenes.lst,gene=CACNA1C

   
// --lookup-gene [pos] chr18:  


// Agreed keywords

   OptionSet * lookup2_opt = par::opt.getOptions("LOOKUP");
   bool useMb = lookup2_opt->isSet("mb") 
     || lookup2_opt->isSet("Mb")
     || lookup2_opt->isSet("MB");

   bool useKb = lookup2_opt->isSet("kb") 
     || lookup2_opt->isSet("Kb")
     || lookup2_opt->isSet("KB");

   if ( useMb && useKb ) 
     error("Cannot specify both Mb and Kb positional queries");

   // A query

   if ( lookup2_opt->isSet("query") )
     {

       NList nl(0);
       NList tlist(0);
       
       vector<string> ids = tlist.deparseStringList( par::lookup2_cmd );
       
       for (int i = 0 ; i < ids.size(); i++)
	 {	   
	   string pquery = ids[i];
	   // Expect format : X=Y
	   size_t p1 = pquery.find("=");
	   string key,val;
	   if (p1==string::npos)
	     {
	       key = pquery;
	       val = "=1";
	     }
	   else
	     {
	       key = pquery.substr(0,p1+1);
	       val = pquery.substr(p1+1);
	   
	       if ( key=="pos=" )
		 {
		   int chr, bp1, bp2;
		   convertPosition(val,chr,bp1,bp2,useKb,useMb);	   
		   val = int2str(chr) + "," + int2str( (int)bp1 ) + "," + int2str( (int)bp2 );
		 }
	     }

	   if ( i == 0 ) 
	     command += key + val;
	   else
	     command += "&"+ key + val;
	 }
     }
   else if ( lookup2_opt->isSet("pos"))
     {
       // Positional query
       // Assume a comma delimited list in format:
       // chr2:737-3993
       
       NList nl(0);
       NList tlist(0);
       vector<string> ids = tlist.deparseStringList( par::lookup2_cmd );

       for (int i = 0 ; i < ids.size(); i++)
	 {
	   // Always requires format X:Y-Z
	   // X should start "chrXX"
	   // Y should be a number
	   // Z should also be a number
	   
	   string pquery = ids[i];	   
	   int chr, bp1, bp2;
	   convertPosition(pquery,chr,bp1,bp2,useKb,useMb);	   
	   pquery = "pos=" + int2str(chr) + "," + int2str( (int)bp1 ) + "," + int2str( (int)bp2 );
	   
	   if ( i == 0 ) 
	     command += pquery;
	   else
	     command += "&"+pquery;
	 }

     }
   else if ( lookup2_opt->isSet("list") || lookup2_opt->isSet("file") )
     {
       checkFileExists(par::lookup2_cmd );
       ifstream IN1( par::lookup2_cmd.c_str() , ios::in );

       command = par::lookup_gene ? "gene=" : "rsid=";
       bool doneFirst = false;
       while ( ! IN1.eof() )
	 {
	   string name;
	   IN1 >> name;
	   if ( name == "" ) 
	     continue;
	   if ( ! doneFirst )
	     {
	       command += name;
	       doneFirst = true;
	     }
	   else
	     command += "," + name;
	 }	 
       IN1.close();
     }
   else if ( lookup2_opt->isSet("qfile") || lookup2_opt->isSet("qlist") )
     {

       checkFileExists(par::lookup2_cmd);
              ifstream IN1( par::lookup2_cmd.c_str() , ios::in );

       bool doneFirst = false;
       while ( ! IN1.eof() )
	 {
	   vector<string> tok = tokenizeLine(IN1);
	   if ( tok.size() == 0 ) 
	     continue;
	   if ( tok.size() != 2 )
	     error("Problem with query file: not 2 columns");
	       
	   string key = tok[0];
	   string val = tok[1];
	   
	   if ( ! doneFirst )
	     {
	       command += key + "=" + val;
	       doneFirst = true;
	     }
	   else
	     command += "&" + key + "=" + val;
	 }	 
       IN1.close();
     }
   else // assume simple gene-name query
     {
       command = par::lookup_gene ? "gene=" : "rsid=";
       
       NList nl(0);
       NList tlist(0);
       vector<string> ids = tlist.deparseStringList( par::lookup2_cmd );
       for (int i = 0 ; i < ids.size(); i++)
	 {
	   if ( ids[i].find("=") != string::npos )
	     error("Badly formed gene name, with equals sign '=' in it -- did you mean to add [query]?");

	   if ( i == 0 ) 
	     command += ids[i];
	   else
	     command += ","+ids[i];
	 }

     }

   
   if ( par::lookup_gene )
     GET_STRING = "GET /plink/genelist.php?";
   else
     GET_STRING = "GET /plink/snplist.php?";
   
   GET_STRING += command;
   //   GET_STRING += "searchtype=gp&chr=6&start=30000000&end=30500000";
   
   cout << "Proposed command = \n";
   cout << GET_STRING << "\n";

   GET_STRING += "\nHTTP/1.0\nContent Length: 10000\nHost: 152.19.78.148\nConnection: close\n\n";
   
   
   cout << "GET_STRING:\n\n" << GET_STRING << "\n";

   
   
   ////////////////////////////////////////////////
   // Make database call
   
   vector<string> tokens = socketConnection( this ,
					     IP_ADDR,
					     PORT_NUM,
					     GET_STRING );
   
   
   ////////////////////////////////////////////////
   // Relay output
   cout << "Output = \n";

   for (int t = 0 ; t < tokens.size() ; t++)
     {
       cout << "token[" << t << "]   =  [ " << tokens[t] << "]\n";
     }
   cout << "\n";
   cout << "-------------------------------\n";
   
#endif
   
}


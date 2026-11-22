

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
#include <sstream>

#include "plink.h"
#include "helper.h"
#include "options.h"
#include "sockets.h"

using namespace std;

#define  PORT_NUM                80     
#define  IP_ADDR    "132.183.161.22"  

void Plink::lookup()
{

#ifdef SKIP
  printLOG("Web-lookup not implemented on this system...\n");
  return;
#else

  printLOG("PLINK-SNP (WGAS SNP annotation courtesy of Patrick Sullivan)\n");

  vector<string> tokens;
  vector< string > gene_list;

    string GET_STRING;

  if ( par::lookup_gene )
    {
      if ( par::lookup_multiple_genes ) 
	{
	  
 	  string query = "";
 	  int c = 0;
	  
 	  string window= int2str(par::lookup_snp_kb_window * 1000);
	  
 	  // Read list, append to query and send
 	  checkFileExists(par::lookup_gene_name);
 	  ifstream S(par::lookup_gene_name.c_str(),ios::in);

	  string cmd;
	  int x = 0;

	  GET_STRING = "GET /~purcell/cgi-bin/gene.pl?win=" + window;

 	  while(!S.eof())
 	    {
 	      string gene;
 	      S >> gene;
 	      if (gene=="")
 		continue;
	      std::string s;
	      std::stringstream out;

	      out << x;
	      s = out.str();
	      cmd += "&gene" + s + "=" + gene;
	      x++;
	      c++;
	      gene_list.push_back(gene);

 	      if ( c>100 ) 
 		error("Please do not send large batch queries to PLINK-SNP");
 	    }
 	  S.close();      
	  GET_STRING += cmd + " HTTP/1.0 \n Content Length: 10000 \nHost: pngu.mgh.harvard.edu\nConnection: close\n\n";
	  tokens = socketConnection( this ,
				     IP_ADDR,
				     PORT_NUM,
				     GET_STRING );

	}
      else
	{
	  string window = int2str( par::lookup_gene_kb_window * 1000 );
	  
	  printLOG("Looking up gene information (and SNPs +/- " 
		   + int2str(par::lookup_gene_kb_window)+" kb)\n");
	  
	  GET_STRING = "GET /~purcell/cgi-bin/gene.pl?win=" + window + "&gene=" + par::lookup_gene_name + " HTTP/1.0 \nHost: pngu.mgh.harvard.edu\nConnection: close\n\n";
	  tokens = socketConnection( this , 
				     IP_ADDR,
				     PORT_NUM,
				     GET_STRING );

	}
    }
  else if ( par::lookup_single_snp )
    {

      string window= int2str(par::lookup_snp_kb_window * 1000);

      printLOG("Looking up SNP information, listing genes within " 
	       + int2str(par::lookup_snp_kb_window)+" kb\n");


      GET_STRING = "GET /~purcell/cgi-bin/snp.pl?win=" + window + "&snp=" + par::lookup_snp + " HTTP/1.0 \nHost: pngu.mgh.harvard.edu\nConnection: close\n\n";
 
      tokens = socketConnection( this , 
				 IP_ADDR,
				 PORT_NUM,
				 GET_STRING );      
    }
  else
    {
      string query = "";
      int c = 0;

      string window= int2str(par::lookup_snp_kb_window * 1000);

      // Read list, append to query and send
      checkFileExists(par::lookup_snp);
      ifstream S(par::lookup_snp.c_str(),ios::in);
      string cmd;
      int x = 0;

      GET_STRING = "GET /~purcell/cgi-bin/snp.pl?win=" + window;
      while(!S.eof())
	{
	  string snp;
	  S >> snp;
	
	  if (snp=="")
	    continue;
	  std::string s;
	  std::stringstream out;

	  out << x;
	  s = out.str();
	  cmd += "&snp" + s + "=" + snp;	
	  x++;
	  c++;
	  if ( c>100 ) 
	    error("Please do not send large batch queries to PLINK-SNP");
	}
      S.close();      
      GET_STRING += cmd + " HTTP/1.0\nHost: pngu.mgh.harvard.edu\nConnection: close\n\n";

      tokens = socketConnection( this ,
      				 IP_ADDR,
      				 PORT_NUM,
      				 GET_STRING );

      printLOG("Looking up SNP information, listing genes within " 
	       + int2str(par::lookup_snp_kb_window)+" kb)\n");

    }
   


  if ( tokens.size() < 25 ) {
      cout << "\n\n";
      if( par::lookup_single_snp ) 
	  cout << par::lookup_snp;
      else
	  cout << par::lookup_gene_name;

      cout << " not found\n";
      cout << "-----------------------------\n";
      return;
  }

   if ( par::lookup_to_file )
    {
      
      ofstream OUT;
      string f = par::output_file_name;
      
      // Need to handle both single and multiple instances

      if ( par::lookup_gene ) 
	{

	  f += ".snp.list";
	  
	  if ( tokens[0] == "-1") 
	    {
	      printLOG("\n\nCould not find gene " 
		       + par::lookup_gene_name 
		       + " in database\n");
	      return;
	    }

	  printLOG("Writing SNP details to [ " + f + " ]\n\n\n");	  
	  OUT.open( f.c_str(), ios::out );

	  bool geneInfo = true;
	  int totalSNPs = 0;
	  int gcount = 0;

	    for (int i=25; i<tokens.size()-1; i++)
	      {
	      if (tokens[i]=="__END")
		{
		  
		  // Skip following line return
		  
		  
		  if ( geneInfo )
		    
		    {
		      // Move from gene information to SNPs
		      
//			  if ( par::lookup_multiple_genes )
//			      OUT << tokens[i+1]; // "\n";
//			      OUT << gene_list[gcount++]; // << "\n";
			  
		      geneInfo = false;
		      continue;
		    }
		  else
		    {
		      // Move from SNPs to gene information
		      // Skip count 
		       ++i;
		      geneInfo = true;
		      printLOG("----------------------------------------------------------\n");
		      printLOG("\n\n");
		      
//		      if ( par::lookup_multiple_genes )
			OUT << "END\n\n";
		      
		      continue;
		    }
		}
	      
	      if ( geneInfo ) 
		{
		  if (tokens[i]=="|n")
		    printLOG("\n");
		  else if ( tokens[i]=="|t")
		    printLOG("\t");
		  else
		    printLOG(tokens[i]+" ");	     
		}
	      else
		{
		  if (tokens[i]=="|n")
		    OUT << "\n";
		  else if ( tokens[i]=="|t")
		    OUT << "\t";
		  else
		    OUT << tokens[i] << " ";	  
		}
	      
	  }
	    OUT.close();
	}
	  else // multiple SNP queries selected
	    {
	      
	      f += ".snp.annot";
	      
	      printLOG("\nWriting SNP details to [ " + f + " ]\n");
	      
	      OUT.open( f.c_str(), ios::out );
	      for (int i=25; i<tokens.size()-1; i++)
		{
		  if (tokens[i]=="|n")
		    OUT << "\n";
		  else if ( tokens[i]=="|t")
		    OUT << "\t";
		  else
		    OUT << tokens[i] << " ";	  
		}      
	      OUT.close();
	    }

    }
    
      else // Standard to LOG (single SNP query)
	{
	  printLOG("\n\n");
	  for (int i=25; i<tokens.size()-1; i++)
	    {
	      if (tokens[i]=="|n")
		printLOG("\n");
	      else if ( tokens[i]=="|t")
		printLOG("\t");
	      else
		printLOG(tokens[i]+" ");
	    }
	}
    
   return;
   
#endif
   
}


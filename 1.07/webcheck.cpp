

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
#include <ctime>

#include "plink.h"
#include "helper.h"
#include "options.h"
#include "sockets.h"

using namespace std;

extern string PVERSION;
extern string PREL;

#define  PORT_NUM                80     
#define  IP_ADDR    "155.52.206.11"
#define  GET_STRING "GET /~purcell/plink/version2.txt HTTP/1.1\nHost: pngu.mgh.harvard.edu\nConnection: close\n\n"


void Plink::webcheck(CArgs & a)
{

#ifdef SKIP
  printLOG("Web-check not implemented on this system...\n");
  return;
#else
  
    
  
  //////////////////////////////////////////
  // First look for a local .pversion file in 
  // the local directory
  
  // Get today's date
  
  time_t curr=time(0);
  string tdstamp = (string)ctime(&curr);      
  string buf; 
  stringstream ss(tdstamp); 
  vector<string> date_tokens; 
  while (ss >> buf)
    date_tokens.push_back(buf);
  string thisDate = date_tokens[0] + date_tokens[1] + date_tokens[2];
  
  bool hasRecord = doesFileExist(".pversion");
  


  ////////////////////////////////////////////////////////
  // Web-based message (but may be cached in local file)

  vector<string> tokens;

  bool connect2web = true;

  printLOG("Web-based version check ( --noweb to skip )\n");


  ////////////////////////////////////////////
  // If we have a record, are we up-to-date?
  
  if ( hasRecord )
    {
      ifstream VER;
      VER.open(".pversion",ios::in);
      string oldDay, oldMonth, oldDate, webVersion;
      VER >> oldDay >> oldMonth >> oldDate;
      
      if ( thisDate == oldDay+oldMonth+oldDate )
	{
	  printLOG("Recent cached web-check found...");
	  connect2web = false;
	  
	  // Read rest of cached web message
	  while ( ! VER.eof() )
	    {
	      string t;
	      VER >> t;
	      if (t=="")
		break;
	      tokens.push_back(t);
	    }	  
	}
      VER.close();
      
    }
  
  
  
  if ( connect2web ) 
    {
      //printLOG("Connecting to web to get version...\n");
      
      tokens = socketConnection( this , 
				 IP_ADDR,
				 PORT_NUM,
				 GET_STRING);      
    }

  bool print = false;
  bool print2 = false;
  bool version_okay = true;

  for (int i=0; i<tokens.size(); i++)
    {

      if (tokens[i]=="END") break;

      if (tokens[i]=="END-MESSAGE")
	{
	  print2=false;
	  continue;
	}

      if (tokens[i]=="WARN")
	{
	  if ( i < tokens.size()-1 ) 
	    {
	      i++;
	      if ( a.find(tokens[i]) )
		{
		  printLOG("\n*** ALERT ***\n*** A warning flag has been set for: "+tokens[i]+
			   "\n*** See http://pngu.mgh.harvard.edu/purcell/plink/warnings.shtml\n");
		  warnings = true;
		}
	    }
	  continue;
	}
      
      
      if (tokens[i]=="FATAL")
	{
	  if ( i < tokens.size()-1 ) 
	    {
	      i++;
	      if ( a.find( tokens[i]) )
		error("A serious warning flag has been set for: "+tokens[i]+
		    "\nPLINK has been instructed to stop"+
 	            "\nPlease see http://pngu.mgh.harvard.edu/purcell/plink/warnings.shtml\n");
	    }
	  continue;
	}


      if (tokens[i]=="MESSAGE-ALL")
	{
	  print2=true;
	  continue;
	}

      // Display any other messages
      // Either conditional on old version (print)
      // or a broadcast to all users (print2)

      if ( ( print && !version_okay) || print2 ) 
	{
	  if (tokens[i]=="\\n")
	    printLOG("\n");
	  else
	    printLOG(tokens[i]+" ");
	}

      // Check version code
      if (tokens[i]=="PLINK-VERSION") 
	{
	  print=true;
	  if ( i < tokens.size() - 1) 
	    {
	      if (tokens[i+1] == PVERSION)
		printLOG(" OK, v"+PVERSION+" is current\n");
	      else
		{
		  printLOG("\n\n          *** UPDATE REQUIRED ***\n\n");
		  printLOG("\tThis version        : "+PVERSION+PREL+"\n");
		  printLOG("\tMost recent version : "+tokens[i+1]+"\n\n");
		  printLOG("Please upgrade your version of PLINK as soon as possible!\n"); 
		  printLOG("  (visit the above website for free download)\n\n");
		  version_okay=false;
		}

	      // Skip the version number
	      i++;
	    }
	}

    }

  // did we get the information we needed?
  if (!print) printLOG("Problem connecting to web\n");

  printLOG("\n");


  ////////////////////////////////////////////////////
  // Create a record that we've checked

  // First line is date-stamp; then simply copy the 
  // whole message
  
  ofstream VER;
  VER.open(".pversion",ios::out);
  VER << date_tokens[0] << " "
      << date_tokens[1] << " "
      << date_tokens[2] << "\n";
  for (int i=0; i<tokens.size(); i++)
    VER << tokens[i] << "\n";
  VER.close();


#endif  
}


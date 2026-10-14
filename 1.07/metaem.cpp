

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2007 Shaun Purcell                  //
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
#include <cmath>
#include <vector>
#include <map>
#include <cassert>

#include "plink.h"
#include "options.h"
#include "helper.h"

#include "genogroup.h"
#include "phase.h"
#include "haplowindow.h"

extern ofstream LOG;
using namespace std;



bool HaploPhase::makeWaplotype(vector<int> & wCounter, vector<int> & wMax)
{
  
  int j = 0;

  while (1)
    {
      if (wCounter[j] < wMax[j])
	{
	  ++wCounter[j];
	  return true;
	}
      else // if this position is a max
	{
	  wCounter[j] = 0;
	  ++j;

	  // Or are we done
	  if (j == actual_nw )
	    return false;
	}
    }

  return true;

}




void HaploPhase::enumeratePhasedWindows(int i)
{
  
  // Form "waplotypes" (haplotypes of windows)
  
  vector<int> wCounter(actual_nw, 0);
  vector<int> wMax(actual_nw, 0);
  
  // Allow at least one move
  wCounter[0] = -1;
  
  // Clear phases for precaution
  hap1[i].clear();
  hap2[i].clear();
  
  for (int w = startWindow; w <= finishWindow ; w++)
    {
      int wc = w - startWindow;
      
      if (windows[w]->hap1[i].size() > 0)
	wMax[wc] = windows[w]->hap1[i].size() - 1;
      else
	{
	  include[i] = false;
	  return;
	}
    }


  //////////////////////////////////////
  // Consider all possible combinations

  while ( makeWaplotype(wCounter, wMax) )
    {

      // But is this a legal window, i.e. given stubs?
      // i.e. scan all intermediate windows and check
      
      bool okay = true;
      
      vector<bool> leftAlignListA;
      vector<bool> leftAlignListB;

      vector<bool> rightAlignListA;
      vector<bool> rightAlignListB;

      for (int w = startWindow; w <= finishWindow; w++)
	{
	  
	  int wc = w - startWindow;
	  
	  // Check left side / left window
	  
	  HaploWindow * currentWindow = windows[w];
	  int currentH1 = currentWindow->hap1[i][ wCounter[wc] ];
	  int currentH2 = currentWindow->hap2[i][ wCounter[wc] ];

	  int l1 = currentWindow->leftStub[ currentH1 ];
	  int l2 = currentWindow->leftStub[ currentH2 ];
	  
	  int r1 = currentWindow->rightStub[ currentH1 ];
	  int r2 = currentWindow->rightStub[ currentH2 ];
	  
	  
	  // Left alignment? 
	  
	  if ( w > startWindow ) 
	    {
	      HaploWindow * leftWindow = windows[w-1];
	      int leftH1 = leftWindow->hap1[i][ wCounter[wc-1] ];
	      int leftH2 = leftWindow->hap2[i][ wCounter[wc-1] ];
	      
	      int ol1 = leftWindow->rightStub[ leftH1 ];
	      int ol2 = leftWindow->rightStub[ leftH2 ];
	      
	      bool leftAlignA = l1 == ol1 && l2 == ol2;
	      bool leftAlignB = l1 == ol2 && l2 == ol1;
	      
	      if ( ! ( leftAlignA || leftAlignB ) )
		{
		  okay = false;
		  break;
		}
	      else
		{
		  leftAlignListA.push_back( leftAlignA );
		  leftAlignListB.push_back( leftAlignB );
		}
	    }
	  else
	    {
	      leftAlignListA.push_back( false );
	      leftAlignListB.push_back( false );
	    }
	  


	  /////////////////////
	  // Right alignment?
	  
	  if ( w < finishWindow ) 
	    {
	      HaploWindow * rightWindow = windows[w+1];
	      int rightH1 = rightWindow->hap1[i][ wCounter[wc+1] ];
	      int rightH2 = rightWindow->hap2[i][ wCounter[wc+1] ];
	      
	      int or1 = rightWindow->leftStub[ rightH1 ];
	      int or2 = rightWindow->leftStub[ rightH2 ];
	      
	      bool rightAlignA = r1 == or1 && r2 == or2;
	      bool rightAlignB = r1 == or2 && r2 == or1;

	      if ( ! ( rightAlignA || rightAlignB ) )
		{
		  okay = false;
		  break;
		}
	      else // save which possible pairings align
		{
		  rightAlignListA.push_back( rightAlignA );
		  rightAlignListB.push_back( rightAlignB );
		}
	    }
	  else
	    {
	      rightAlignListA.push_back( false );
	      rightAlignListB.push_back( false );
	    }
	}
      
      

      /////////////////////////////////////////////////////////
      // For legal combinations, enumerate possible waplotypes

      if ( okay )
	{
	  
	  vector<bool> flipable(actual_nw,false);
	  int nflip = 0;	  

	  bool firstSkipped = false;
	  
	  // Consider all A/B and B/A pairings of haplotypes for 
	  // windows other than the first (i.e. we will be lining up 
	  // haplotypes in hap1/hap2 meaningfully then)

//	  if ( startWindow > 0 ) firstSkipped = true;
	  
	  for (int w = startWindow; w <= finishWindow; w++)
	    {
	      int wc = w - startWindow;

	      HaploWindow * currentWindow = windows[w];
	      int currentH1 = currentWindow->hap1[i][ wCounter[wc] ];
	      int currentH2 = currentWindow->hap2[i][ wCounter[wc] ];
	      
	      if ( currentH1 != currentH2 )
		{
		  if ( firstSkipped )
		    {
		      flipable[wc] = true;
		      nflip++;
		    }
  //		  firstSkipped = true;
		}		 
	    }
	  
	  
	  int nperm = (int)pow(2.0,(double)nflip);
	  
	  unsigned int h=0;
	  
	  while (h<nperm)
	    {

	      vector<bool> f1;
	      
	      unsigned int p=1;
	      for (int flip=0; flip<nflip; flip++)
		{
		  if ( h & p )
		    f1.push_back(false);
		  else
		    f1.push_back(true);
		  p <<= 1;
		}
	      
	      vector<int> waplotype1;
	      vector<int> waplotype2;

	      int f2 = 0;
	      for (int w = startWindow; w <= finishWindow ; w++)
		{

		  int wc = w - startWindow;

		  HaploWindow * currentWindow = windows[w];

		  int currentH1 = currentWindow->hap1[i][ wCounter[wc] ];
		  int currentH2 = currentWindow->hap2[i][ wCounter[wc] ];
		  
		  if ( flipable[wc] )
		    {
		      if ( f1[f2++] )
			{
			  int tmp = currentH1;
			  currentH1 = currentH2;
			  currentH2 = tmp;
			}
		    }
		  

		  waplotype1.push_back( currentH1 );
		  waplotype2.push_back( currentH2 );

		}
	      
	      
	      ///////////////////////////////
	      // Test for ultimate alignment

	      bool allAligned = true;

	      // Skip first window
	      for ( int w = startWindow+1; w <= finishWindow; w++)
		{
		  int wc = w - startWindow;

		  HaploWindow * currentWindow = windows[w];
		  
		  int l1 = currentWindow->leftStub[ waplotype1[wc] ];
		  int l2 = currentWindow->leftStub[ waplotype2[wc] ];
		  
		  HaploWindow * leftWindow = windows[w-1];

		  int ol1 = leftWindow->rightStub[ waplotype1[wc-1] ];
		  int ol2 = leftWindow->rightStub[ waplotype2[wc-1] ];
		  
		  bool leftAlign = l1 == ol1 && l2 == ol2;
		  
		  if ( ! leftAlign ) 
		    allAligned = false;
		}
	      

	      if ( allAligned )
		{


		  // Have we already seen these two waplotypes?
		  
		  map< vector<int>, int>::iterator wi = hapmap.find( waplotype1 );
		  
		  // First waplotype 


		  if ( wi == hapmap.end() )
		    {

		      hapmap.insert( make_pair( waplotype1 , nh ));

		      vector<bool> thisHaplotype;
		      for(int w= startWindow; w <= finishWindow; w++)
			{
			  int wc = w - startWindow;

			  HaploWindow * currentWindow = windows[w];
			  int start = par::haplo_plem_overlap;
			  if (wc==0) start = 0;
			  for (int s = start; s < currentWindow->ns; s++)
			    thisHaplotype.push_back( currentWindow->hap[waplotype1[wc]][s]);
			}
		      
		      hap.push_back( thisHaplotype );
		      hapmapb.insert( make_pair( thisHaplotype, nh ));
		      hapi.push_back( waplotype1 );

		      hap1[i].push_back( nh );

		      nh++;
		    }
		  else
		    {
		      hap1[i].push_back( wi->second );
		      hapmapb.insert( make_pair( hap[wi->second], wi->second ));
		    }

	
		  wi = hapmap.find( waplotype2 );
		  
		  if ( wi == hapmap.end() )
		    {
		      hapmap.insert( make_pair( waplotype2 , nh ));
		      
		      vector<bool> thisHaplotype;
		      for(int w=startWindow; w <= finishWindow; w++)
			{
			  int wc = w - startWindow;

			  HaploWindow * currentWindow = windows[w];
			  int start = par::haplo_plem_overlap;
			  if (wc==0) start = 0;
			  for (int s = start; s < currentWindow->ns; s++)
			    thisHaplotype.push_back( currentWindow->hap[waplotype2[wc]][s]);
			}
		      hap.push_back( thisHaplotype );
		      hapmapb.insert( make_pair( thisHaplotype, nh ));
		      hapi.push_back( waplotype2 );
		      hap2[i].push_back( nh );
		      
		      nh++;
		    }
		  else
		    {
		      hap2[i].push_back( wi->second );
		      hapmapb.insert( make_pair( hap[wi->second], wi->second ));
		    }

		}

			  
	      // Consider next waplotype
	      h++;
	      
	    }
       	
	} // end of 'add if legal' condition      

    } // consider next possible legal set
  
  
  // If more than a single phased set exists declare ambigious

  ambig[i] = hap1[i].size() > 1;
  include[i] = hap1[i].size() > 0;
  
  if ( par::haplo_plem_follow && i == par::haplo_plem_follow_ind )
    {
      VPHASE << "Added " << hap1[i].size() 
	     << " phases for followed individual\n";      
    }

}





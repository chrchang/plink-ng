


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

#include "plink.h"
#include "options.h"
#include <cmath>

bool isAncestorOf(Individual *indx, Individual * f);
int mCount(Individual * indx, Individual *f);

void listAllAncestors(Individual * a, set<Individual*> & anclist, int d)
{
  
  anclist.insert( a );

  if ( a->pm == NULL && a->pp == NULL )
    return;
  
  if ( a->pp )
    {
      anclist.insert( a->pp );
      listAllAncestors( a->pp , anclist , d+1 );
    }
  
  if ( a->pm )
    {
      anclist.insert( a->pm );
      listAllAncestors( a->pm , anclist , d+1 );
    }
  
  return;
  
}


double genrel(Individual * a, Individual * b) 
{

  //  cout << "GR for " << a->iid << " and " << b->iid << "\n";

  double g = 0;

  // Same person?

  if ( a == b )
    return 1;
	  
  // Are both individuals founders or in different families?
  
  if ( a->fid != b->fid )
    return 0;

  if ( a->founder && b->founder )
    return 0;
  
  
  // Assuming no inbreeding, find the nearest common ancestor
  
  // For each possible individual, store the number of meioses that 
  // separate a from k and b from k (in an int2, -1 for not a common 
  // ancestor)

  map<Individual*,int2> nca;
  
  // Start with ancestors of A
  
  set<Individual*> ancestorsA;
  set<Individual*> ancestorsB;

  listAllAncestors(a,ancestorsA,0);
  listAllAncestors(b,ancestorsB,0);
  
  multiset<Individual*> commonAncestors;
  set<Individual*>::iterator i = ancestorsA.begin();
  while( i != ancestorsA.end() )
    {
      commonAncestors.insert( *i );
      ++i;
    }
  i = ancestorsB.begin();
  while( i != ancestorsB.end() )
    {
      commonAncestors.insert( *i );
      ++i;
    }

  // Any individuals represented twice?
  set<Individual*> mrca;
  multiset<Individual*>::iterator j = commonAncestors.begin();
  while ( j != commonAncestors.end() )
    {
      if ( commonAncestors.count( *j ) == 2 )
	mrca.insert( *j );
      ++j;
    }

//    cout << "sizes = " << ancestorsA.size() << " " << ancestorsB.size() << "\n";
//    cout << "MRCA for " << a->fid << " " << a->iid << " / " << b->iid << "\n";

//    i = mrca.begin();
//    while( i != mrca.end() )
//      {
//        cout << (*i)->fid << " " << (*i)->iid << "\n";
//        ++i;
//      }
//    cout << "\n";

  
  // Iterate through common ancestors, finding # of 
  // meioses back to the two founders individuals
  
  i = mrca.begin();
  while( i != mrca.end() )
    {      
      int2 m( (*i)->countMeioses(a) , (*i)->countMeioses(b) );
      nca.insert( make_pair( *i, m ) );      
      ++i;
    }
  
  int kmin = 9999;
  i = mrca.begin();
  while( i != mrca.end() )
    {      
      
      int2 m( (*i)->countMeioses(a) , (*i)->countMeioses(b) );      
      nca.insert( make_pair( *i, m ) );      
      
      if ( m.p1 + m.p2 < kmin ) 
	kmin = m.p1 + m.p2;
      
      ++i;
    }


  //////////////////////
  // Calculate 'g'

  map<Individual*,int2>::iterator k = nca.begin();
  while( k != nca.end() )
    {      
      int2 m = k->second;
      int m2 = m.p1 + m.p2;
      
      if ( m2 == kmin )
	{
	  g += pow(0.5,m2);
	  //cout << "adding " <<  k->first->iid << " : " << m.p1 << " + " << m.p2 << "\n";
	}
      ++k;
    }

  return g;

}




int Individual::countMeioses(Individual *f)
{
  
  if ( isAncestorOf(this,f ) ) 
    return mCount(this,f);
  else if ( isAncestorOf( f,this ) ) 
    return mCount(f,this);
  else
    return 0;      
}



int mCount(Individual * indx, Individual *f)
{
  
  vector<Individual*> inds;
  vector<bool> checked;
  bool finished = false;
  int nm = 0;
  
  
  // Add self to list
  
  inds.push_back(indx);
  checked.push_back(false);
  
  while (!finished)
    {
      
      // Check list for a match
      // needs changing if inbreeding
      
      for (int i = 0 ; i < inds.size() ; i++)
	{ 
	  if (inds[i] == f) return nm; 
	}
      
      // Increment meioses counter
      nm++;
      
       // Add children of unchecked inds

      int already = inds.size();

      for (int i = 0 ; i < already ; i++)
	{
	  if (!checked[i])
	    {
	      for (int j = 0 ; j < inds[i]->kids.size() ; j++)
		{

		  inds.push_back(inds[i]->kids[j]);
		  checked.push_back(false);
		  		  
		}
	      checked[i] = true;
	    }
         }


       // All done?

       finished = true;

       for (int i = 0 ; i < inds.size() ; i++)
         if ( checked[i] == false ) finished = false;

       // loop back
     }


  return nm;
}


bool isAncestorOf(Individual *indx, Individual * f)
{

  vector<Individual*> inds;
  vector<bool> checked;
  bool finished = false;
  int nm = 0;

  // Add self to list
  inds.push_back(indx);
  checked.push_back(false);

  while (!finished)
    {
      // Check list for a match
      for (int i = 0 ; i < inds.size() ; i++)
          if (inds[i] == f) return true;


      // Increment meioses counter
      nm++;

      // Add children of unchecked inds
      int already = inds.size();
      for (int i = 0 ; i < already ; i++)
        {
          if (!checked[i])
            {
              for (int j = 0 ; j < inds[i]->kids.size() ; j++)
                {
                  inds.push_back( inds[i]->kids[j] );
                  checked.push_back(false);
                }
              checked[i] = true;
            }
        }


      // All done?
      finished = true;
      for (int i = 0 ; i < inds.size() ; i++)
        if ( checked[i] == false ) finished = false;
      // loop back
    }
  return false;
}

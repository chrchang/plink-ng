

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <algorithm>

#include "nlist.h"

vector<string> NList::deparseStringList(string input)
{
  return tokenize(input);
}

vector<int> NList::deparseNumberList(string input)
{

  // We have a string that should contain integers, possibily
  // separated by "," and "-" characters

  // Use 1-N coding up until the last moment, instead of 
  // standard 0..N-1 coding

  firstWord = "1";
  lastWord = int2str(maxcat);

  vector<string> tok = tokenize(input);
  vector<int> nums;
  for (int i=0; i<tok.size(); i++)
    {
      int j;
      if ( tok[i] == range_string )
	nums.push_back(-1);
      else
	{
	  if ( ! from_string<int>( j, tok[i] , std::dec))
	    continue;
	  nums.push_back(j);
	}
    }
  return expandNumberList(nums);
}


vector<int> NList::deparseStringList(string input, map<string,int> * mapping)
{
  
  // Assume map contains a 0..N-1 coding in the map
  
  map<string,int>::iterator im = mapping->begin();
  while ( im != mapping->end() )
    {
      int i = im->second;
      
      if ( i == 0 )
	firstWord = im->first;
      
      if ( i == maxcat - 1 ) 
	lastWord = im->first;

      ++im;
    }


  ////////////////////////////////////////////////////////////////
  // Convert string codes to numbers, then call deparseNumberList
  
  // But incrememnt to 1..N coding here, i.e. so that it is 
  // similar process to the human-friendly 1..N coding for a
  // numeric list

  vector<string> tok = tokenize(input);
  vector<int> nums;
  for (int i=0; i<tok.size(); i++)
    {
      map<string,int>::iterator im = mapping->find(tok[i]);
      if ( im != mapping->end() )
	nums.push_back(im->second + 1);
      else if ( tok[i] == range_string )
	nums.push_back(-1);
      else 
	error("Cannot find value: " + tok[i] + "\n");
    }
  return expandNumberList(nums);  
}



vector<int> NList::expandNumberList(vector<int> & nlist)
{

  // Convert 
  vector<int> n;
  
  //////////////////
  // Check ranges

  for (int i=0; i<nlist.size(); i++)
    {
      if ( nlist[i] == -1 ) 
	{
	  if ( nlist[i-1] > nlist[i+1] )
	    {
	      int tmp = nlist[i+1];
	      nlist[i+1] = nlist[i-1];
	      nlist[i-1] = tmp;
	    }
	}
    }
  

  //////////////////
  // Expand ranges
  
  for (int i=0; i<nlist.size(); i++)
    {
      
      if ( nlist[i]>0 ) 
 	{
	  // Only add valid codes
	  if ( nlist[i] <= maxcat )
	    n.push_back(nlist[i]);	  
 	}
      else
	{
	  int start = nlist[i-1]+1;
 	  int end = nlist[i+1]-1;
	  if ( end > maxcat )
	    end = maxcat;

 	  for (int j=start; j<=end; j++)
 	    n.push_back(j);
 	}
    }

  // Sort and uniquify
  stable_sort(n.begin(),n.end());
  vector<int>::iterator ne = unique(n.begin(),n.end());
  n.erase(ne, n.end());

  // Shift to 0..N-1 coding
  for (int i=0; i<n.size(); i++)
    --n[i];


  // Or is it a negative complement?
  if ( negmode )
    {

      vector<bool> k(maxcat,false);
      for (int i=0; i<n.size(); i++)
	k[ n[i] ] = true;
      
      vector<int> n2 = n;
      n.clear();
      
      for (int i=0; i<maxcat; i++)
	{
	  if ( ! k[i] ) 
	    {
	      n.push_back(i);
	    }
	}
    }
 
  return n;

}


vector<string> NList::tokenize(string s)
{

  vector<string> t;

  string word = "";
  bool range_set = false;
  bool number_given = false;
  
  for (int i=0; i<s.length(); i++)
    {      

      // Skip spaces
      if (s[i] == ' ') 
	continue;
      
      if ( s[i] == delimit_char || s[i] == range_char )
	{

	  // If we've built up a meaningful word,
	  // we can add it now

	  if ( word != "" )
	    {
	      t.push_back(word);
	      word = "";
	      number_given = true;
	      range_set = false;
	    }
	  
	  
	  // Start a range? 
	  if ( s[i]== range_char )
	    {

	      // Can't do twice
	      if (range_set)
		error("Invalid list: " + s + "\n");
	      
	      range_set = true;
	      
	      // But was a start point given? Make one 
	      // if not

	      if ( ! number_given )
		{
		  t.push_back("1");
		  word="";
		}
	      
	      number_given = false;
	      
	      t.push_back(range_string);
	      
	    }
	  
	  // And if the "-" is also the last 
	  // character, we need an end point
	  
	  if ( s[i]== range_char && i == s.length()-1 )
	    t.push_back( lastWord );
	  
	  // Here, fill in here
	  if ( s[i]== delimit_char )
	    {
	      number_given = false;
	      if ( range_set )
		{
		  range_set = false;
		  t.push_back( lastWord );
		}	  
	    } 
	}
      else if ( i == s.length()-1 )
	{
	  // Is this a last number?
	  t.push_back(word+s[i]);
	  number_given = true;
	}
      else
	word += s[i];
    }

  return t;
}



//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2008 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#ifndef __NLIST_H__
#define __NLIST_H__

#include <string>
#include <vector>
#include <iostream>

#include "helper.h"
#include "options.h"

using namespace std;
class NList{
    
  vector<string> tokenize(string);
  vector<int> expandNumberList(vector<int> &);
  int maxcat;
  bool negmode; // exclude list from 1..maxcat

  string firstWord;
  string lastWord;

  char range_char;
  string range_string;
  
  char delimit_char;
  
 public:

  NList(int n, bool nmode = true )
    { 
      range_char = par::range_delimiter[0];
      range_string = par::range_delimiter;
      delimit_char = ',';
      maxcat = n; 
      negmode = ! nmode;
    }
  void setRangeChar(string s)
    {
      range_char = s[0];
      range_string = s;
    }

  void setDelimiter(string s)
    {
      delimit_char = s[0];
    }


  vector<int> deparseNumberList(string);
  vector<int> deparseStringList(string, map<string,int> *);
  vector<string> deparseStringList(string);
};

#endif

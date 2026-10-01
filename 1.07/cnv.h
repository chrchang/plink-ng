

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#ifndef __CNV_H__
#define __CNV_H__

class CNVIndivReport {
 public:

  CNVIndivReport() 
    {
      n = count = count_baseline = 0;
      segCount = 0;
      t1 = t2 = t3 = t4 = t5 = t6 = t7 = t8 = 0;
    }
  
  void calculateResults();

  // Total number of individuals looked at
  int n;

  // Number of individuals with an event
  int count;

  // Number of individuals with a baseline-region event
  int count_baseline;

  // Number of segments in this group
  int segCount;
  
  // Test results
  // 1) # segs; 2) #inds w 1+; 3) total kb; 4) mean seg kb, 5) gene-count, 
  // 6) gene-count-1+, 7) base-line gene-count
  
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  double t7;
  double t8;

  // if needed...
  double t9;
  double t10;
  double t11;
  double t12;

};


// Helper functions

bool intersects(set<Range>&,set<Range>&,int,int,int);
int count_intersects(set<Range>&,int,int,int);
double weighted_count_intersects(set<Range>&,int,int,int);
vector<int> segmentCountCaseControls(Plink*,bool);
set<Segment> allSegmentsIntersecting(Range &);


#endif

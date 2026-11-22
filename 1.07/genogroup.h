

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2006 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#ifndef __GENOGROUP_H_
#define __GENOGROUP_H__


class MultiLocusGenotype {
 public:

  vector<bool> g;
  int count;
  int reference;
  vector<bool> skip;

  bool operator<(const MultiLocusGenotype & b) const {
    for (int i=0; i<g.size(); i++)
      if (g[i] != b.g[i])
	return g[i];
    return false;
  }

  bool operator==(const MultiLocusGenotype & b) const {
    for (int i=0; i<g.size(); i++)
      if (g[i] != b.g[i])
	return false;
    return true;
  }

};
namespace std {
  template<> class less<MultiLocusGenotype*> {
  public:
    bool operator()(MultiLocusGenotype const* p1,
		    MultiLocusGenotype const* p2) const {
      if (!p1)
	return true;
      if (!p2)
	return false;

      if (p1->g < p2->g)
	return true;

      return false;
    }
  };
};


#endif

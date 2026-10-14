

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2006 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#include <iostream>

#include "plink.h"
#include "perm.h"
#include "options.h"
#include "helper.h"
#include "crandom.h"



void Perm::preGeneDrop()
{

  // Note -- minor issue, this routine ignores 
  // issue of linkage between sibs 
  
  // Idea: to use standard case/control or QT test, but permute only
  // transmissions from founders, and to all offspring: to give a
  // within-family test.

  
  // Set up parent-offspring structure, but not related to 
  // nuclear families



  // par::perm_genedrop         
  // If true, perform gene-dropping permutation instead of label-swapping

  // par::perm_genedrop_founders
  // If true, non-founder parents always drop one of their gene-dropped alleles
  // If false, non-founder parents always drop one of their true alleles

  // par::perm_genedrop_parents
  // If true, we also perform a label-swapping permutation within all parents

  // par::perm_genedrop_sibships
  // If true, we also perform a label-swapping permutation within all full 
  // sibships without parents


  map<string,Individual*> fnd;
  
  // Link up parents and offspring
  P.linkRelateds(idmap, fnd);
  
  P.printLOG("Allocated family structure for gene-dropping\n");
  
  
  // Set initial permutation structure -- no label-swapping
  for (int i=0; i<P.sample.size(); i++)
    P.sample[i]->sol = -1;
  

  // Label-swapping cluster count
  int cc=0;
  int cc_par=0;

  map<string,int> parent;
  map<string,int> parent_pat;
  map<string,int> parent_mat;
  
  // Set up clusters for within-parent label-swapping permutation
  if (par::perm_genedrop_parents)
    {
      
      // Parents must be pairs for simple nuclear families
      // i.e. watch out for half-sib relations
      
      for (int i=0; i<P.n; i++)
	{

	  Individual * person = P.sample[i];
	  
	  // Only consider non-founders
	  if ( fnd.find(person->fid+"_"+person->iid) == fnd.end() )
	    {

	      string spat = person->pp->fid+"_"+person->pp->iid;
	      string smat = person->pm->fid+"_"+person->pm->iid;
	      string pair = spat+" x "+smat;
	      
	      // If this parent pair has not previously featured, 
	      // AND if neither parent has previously featured in 
	      // any other pairing, then make this parental set 
	      // a pair
	      
	      if ( parent.find(pair) == parent.end() && 
		   parent_pat.find(spat) == parent_pat.end() && 
		   parent_mat.find(smat) == parent_mat.end() )
		{
		  person->pp->sol = cc;
		  person->pm->sol = cc;
		  cc++;
		  
		  parent.insert(make_pair(pair,cc));
		  parent_pat.insert(make_pair(spat,cc));
		  parent_mat.insert(make_pair(smat,cc));
		}
	    } 
	}
      P.printLOG("Allocated "+int2str(cc)+" clusters for within-parent permutation\n");  
      cc_par = cc;
    }
  
  // Set up clusters for within-sibship permutation
  if (par::perm_genedrop_sibships)
    {

      // i.e for individuals for whom pat and mat != 0 but the parent is 
      // no longer in the dataset (i.e. removed for low genotyping, as it
      // was a dummy parent).
      
      map<string,int> sibs;
      
      for (int i=0; i<P.n; i++)
	{

	  Individual * person = P.sample[i];
	  
	  // Only consider true non-founders i.e. unlike the above
	  // section, which considers only parents, but only those 
	  // for whom we do not have 2 parents
	  
	  if ( ! person->founder && 
	       fnd.find(person->fid+"_"+person->iid) != fnd.end() &&
	       parent_pat.find(person->fid+"_"+person->iid) == parent_pat.end() &&
	       parent_mat.find(person->fid+"_"+person->iid) == parent_mat.end() ) 
	    {
	      
	      string pair = person->fid+"_"+person->pat+"_"+person->mat;
	      
	      // If we haven't seen this sibship before, add a new cluster code
	      map<string,int>::iterator sit = sibs.find(pair);
	      if (sit == sibs.end())
		{
		  person->sol = cc;
		  sibs.insert(make_pair(pair,cc));
		  cc++;
		}
	      else
		{
		  // ...otherwise, assign to existing one
		  person->sol = sit->second;
		}
	      
	    }
	}

      P.printLOG("Allocated "+int2str(cc-cc_par)+" clusters for within-sibship permutation\n");  
    }

  
  // Label-swapping permutation of all unrelated individuals?
  // everybody else, who is a family size 1 
  // this means 
  if (par::perm_genedrop_unrel)
    {
      map<string,int> unrel;
      for (int i=0; i<P.n; i++)
	{
	  string f = P.sample[i]->fid;
	  
	  if ( unrel.find(f) == unrel.end() )
	    {
	      unrel.insert(make_pair(f,1));
	    }
	  else
	    {
	      (unrel.find(f)->second)++;
	    }
	}
      
      for (int i=0; i<P.n; i++)
	cout << P.sample[i]->fid << "\t" << unrel.find(P.sample[i]->fid)->second << "\n";
      
      
      // If no parents, assign unique cluster
//       if (fnd.find(person->fid+"_"+person->iid) != fnd.end())
// 	person->sol = cc;

      
      P.printLOG("Allocated cluster for between-founder permutation\n");  
      
    }
}


void Perm::geneDrop()
{

  // Transmissions
  vector<bool> pat(P.n,true);
  vector<bool> mat(P.n,false);
  vector<bool> done(P.n,false);

  // Consider each individual
  for (int i=0; i<P.n; i++)
    {
      Individual * person = P.sample[i];

      // 1. For non-founders, permute which alleles they inherited
      if ( !person->founder ) 
	{
	  pat[i] = mat[i] = false;
	  if (CRandom::rand() > 0.5) pat[i] = true;
	  if (CRandom::rand() > 0.5) mat[i] = true;
	}
    }
  
  
  // Now we have constructed the gene-dropping matrix, 
  // we proceed to consider each SNP at a time

  for (int l=0; l<P.nl_all; l++)
    {
      fill(done.begin(), done.end(), false);
      for (int i=0; i<P.n; i++)
	if (P.sample[i]->founder) 
	  dropAlleles(P,P.sample[i],i,l,pat,mat,done,idmap);
    }

}


void Perm::dropAlleles(Plink & P,
		       Individual * person, 
		       int i,
		       int l,
		       vector<bool> & pat, 
		       vector<bool> & mat,
		       vector<bool> & done,
		       map<Individual*,int> & idmap)
{
    
  // If founder, leave genotype as is; also, if either parent has
  // missing genotype data, then do not permute
  
  vector<bool>::iterator s1;
  vector<bool>::iterator s2;

  bool pat1, pat2;
  bool mat1, mat2;

  if (par::SNP_major)
    {
      s1 = P.SNP[l]->one.begin()+i;
      s2 = P.SNP[l]->two.begin()+i;
      
      pat1 = P.SNP[l]->one[person->ip];
      pat2 = P.SNP[l]->two[person->ip];

      mat1 = P.SNP[l]->one[person->im];
      mat2 = P.SNP[l]->two[person->im];           
    }
  else
    {
      s1 = person->one.begin()+l;
      s2 = person->two.begin()+l;
      
      pat1 = P.sample[person->ip]->one[l];
      pat2 = P.sample[person->ip]->two[l];

      mat1 = P.sample[person->im]->one[l];
      mat2 = P.sample[person->im]->two[l];
    }

  
  if ( ! ( person->founder ||       // founder
	   ( pat1 && !pat2 ) ||     // pat missing
	   ( mat1 && !mat2 ) ||     // mat missing
	   ( (*s1)   && ! *s2   ) ) )    // self missing
    
    {
      
      // For pat/mat :  
      //   false = paternal/slot1, 
      //   true = maternal/slot2

      // i.e. if parent is heterozygous, then pat T/F says which
      // allele to take (F/T), otherwise just take the homozygous
      // allele

      bool d1 = false;
      bool d2 = false;
      
      // Is father heterozygous?
      if ( pat1 != pat2 )
	d1 = pat[i];
      else
	d1 = pat1;
      
      // Is mother heterozygous?
      if ( mat1 != mat2 )
	d2 = mat[i];
      else
	d2 = mat1;
      
      
      // Set new genotype: FF, FT or TT?
      // (Missing will be left as is)

      
      if ( (!d1) && (!d2) ) 
	{
	  *s1 = false;
	  *s2 = false;
	}
      else if ( d1 != d2 )
	{
	  *s1 = false;
	  *s2 = true;
	}
      else if ( d1 && d2 ) 
	{
	  *s1 = true;
	  *s2 = true;
	}
    
          
      done[i]=true;
      
    }

    
  // Now also update any kids of this person that still need doing
  for (int k=0; k<person->kids.size(); k++)
    if (!done[idmap.find(person->kids[k])->second]) 
      dropAlleles(P,person->kids[k],person->ikids[k],l,pat,mat,done,idmap);
  
  return;
}





void Plink::linkRelateds(map<Individual*,int> & idmap,
			 map<string,Individual*> & fnd)
{
  
  map<string,Individual*> imap; 
  map<Individual*,int> imap2;
  map<string,Individual*>::iterator iit;
  map<Individual*,int>::iterator iit2;
  
  // Populate map, clear any existing family-related information
  for (int i=0; i<n; i++)
    {
      imap.insert(make_pair(sample[i]->fid+"_"+
			    sample[i]->iid,
			    sample[i]));
      imap2.insert(make_pair(sample[i],i));
      idmap.insert(make_pair(sample[i],i));
      sample[i]->kids.clear();
      sample[i]->ikids.clear();
      sample[i]->family = NULL;
    }
  
  
  // Link up parents and offspring

  for (int i=0; i<n; i++)
    {
      Individual * person = sample[i];
      
      if (person->founder)
	{
	  person->pp = person->pm = NULL;
	  person->ip = person->im = -1;
	  fnd.insert(make_pair(sample[i]->fid+"_"+sample[i]->iid,sample[i]));
	}
      else
	{
	  // Father (if does not exist, treat as founder)
	  iit = imap.find(person->fid+"_"+person->pat);
	  if (iit == imap.end())
	    {
	      person->pp = NULL;
	      fnd.insert(make_pair(sample[i]->fid+"_"+sample[i]->iid,sample[i]));
	    }
	  else
	    {
	      person->pp = iit->second;
	      iit2 = imap2.find(iit->second);
	      person->ip = iit2->second;
	    }

	  // Mother (if does not exist, treat as founder)
	  iit = imap.find(person->fid+"_"+person->mat);
	  if (iit == imap.end())
	    {
	      person->pm = NULL;
	      fnd.insert(make_pair(sample[i]->fid+"_"+sample[i]->iid,sample[i]));
	    }
	  else
	    {
	      person->pm = iit->second;	  
	      iit2 = imap2.find(iit->second);
	      person->im = iit2->second;
	    }
	  
	  // Otherwise, add this person as a child of mother and father
	  if ( ! (person->pp == NULL || person->pm == NULL ) )
	    {
	      person->pp->kids.push_back(person);
	      person->pm->kids.push_back(person);
	      
	      person->pp->ikids.push_back(i);
	      person->pm->ikids.push_back(i);
	      
	    }
	}
  
  }
}

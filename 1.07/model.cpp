

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#include "model.h"
#include "options.h"
#include "helper.h"
#include "phase.h"
#include "stats.h"

#include <cmath>

Model::Model()
{    
  np = nind = 0;  
  haploid.resize(0);
  xchr.resize(0);
  order.clear();
  sex_effect = false;
  all_valid = true;
  has_snps = true;
  testParameter = 1; // Permutation test parameter
  
  // Automatically add intercept now
  label.push_back("M"); // Intercept
  type.push_back( INTERCEPT );
  order.push_back(0);


  ///////////////////////////
  // Default additive coding

  mAA = 0;
  mAB = 1;
  mBB = 2;


  ///////////////////////////
  // Set X chromosome coding

  if ( par::xchr_model == 1 )
    {
      mA = 0; 
      mB = 1;
    }
  else if ( par::xchr_model == 2 )
    {
      mA = 0;
      mB = 2;
    }
  else if ( par::xchr_model > 2 )
    {
      mA = 0;
      mB = 1;
    }

}

void Model::setDominant()
{
  mAA = 0;
  mAB = 1;
  mBB = 1;
  mA = 0;
  mB = 1;
}

void Model::setRecessive()
{
  mAA = 0;
  mAB = 0;
  mBB = 1;

  // No haploid effect
  mA = mB = 0;
}

void Model::addSexEffect()
{
  sex_effect = true;
  type.push_back( SEX );
  order.push_back(0);
}

bool Model::isSexInModel()
{
  return sex_effect;
}

void Model::hasSNPs(bool b)
{
  has_snps = b;
}

void Model::setMissing()
{
  
  // Fill in missing data with existing pattern 
  // and also optional per-test missingness

  miss.clear();
  miss.resize(P->n,false);
  for (int i=0; i<P->n; i++)
    if ( P->sample[i]->missing || 
	 P->sample[i]->missing2 ) miss[i] = true;
  
}

vector<bool> Model::getMissing()
{
  return miss;
}

void Model::yokeMissing(Model * m)
{
  //
}

void Model::setMissing(vector<bool> & include)
{
  
  // Fill in missing data with existing pattern 
  if ( include.size() != P->n )
    error("A problem in setMissing()\n");

  miss.resize(P->n,false);
  for (int i=0; i<P->n; i++)
    if ( P->sample[i]->missing || ! include[i] ) miss[i] = true;
  
}


void Model::addAdditiveSNP(int a)
{

  if ( ! has_snps )
    error("Cannot add SNP to this MODEL");

  additive.push_back(a);

  if ( par::chr_sex[P->locus[a]->chr] )
    xchr.push_back(true);
  else
    xchr.push_back(false);

  if ( par::chr_haploid[P->locus[a]->chr] )
    haploid.push_back(true);
  else
    haploid.push_back(false);

  type.push_back( ADDITIVE );
  order.push_back( additive.size() - 1 );

}

void Model::addDominanceSNP(int d)
{
  if ( ! has_snps )
    error("Cannot add SNP to thie MODEL");

  dominance.push_back(d);

  type.push_back( DOMDEV );
  order.push_back( dominance.size() - 1 );

}

void Model::addCovariate(int c)
{
  covariate.push_back(c);      

  type.push_back( COVARIATE );
  order.push_back( covariate.size() - 1 );

}

void Model::addHaplotypeDosage(set<int> & h)
{
  haplotype.push_back(h);

  type.push_back( HAPLOTYPE );
  order.push_back( haplotype.size() - 1 );

}

void Model::addInteraction(int a, int b)
{
  int2 i;
  i.p1 = a;
  i.p2 = b;
  interaction.push_back(i);

  type.push_back( INTERACTION );
  order.push_back( interaction.size() - 1 );

}

void Model::buildDesignMatrix()
{
  
  // Build X matrix (including intercept)
  // Iterate a person at a time, entering only 
  // valid rows into X (non-missing); also build Y
  // at the same time

  if ( has_snps && par::SNP_major )
    error("Internal error: must be individual-major to perform this...\n");

    
  ///////////////////////
  // Number of parameters
  
  // Standard variables
  // Note: 'additive' really means 'main effect' here, 
  //       i.e. which can also be coded recessive or dominant
  //       i.e. the distinction is between the 2df general model
  
  np = 1 
    + additive.size() 
    + dominance.size() 
    + haplotype.size()
    + covariate.size() 
    + interaction.size(); 
  
  // Sex effect?
  if ( sex_effect )
    np++;
  
  // QFAM variables
  //  if (par::QFAM_total) np++;
  //  else if (par::QFAM_between || par::QFAM_within1 || par::QFAM_within2) np +=2;
  
  if (par::QFAM_total
      || par::QFAM_between 
      || par::QFAM_within1 
      || par::QFAM_within2) 
    {
      np++;
      
      type.push_back( QFAM );
      order.push_back( 0 );
    }


  
  
  ///////////////////////////
  // Consider each individual

  for (int i=0; i < P->n; i++)
    {

      Individual * person = P->sample[i];

      // Ignore if missing phenotype, or the user set this to missing
      if ( miss[i] ) continue;
    
      /////////////////////////////
      // 0) Intercept      
      // 1) Main effects of SNPs 
      // 2) Dominance effects of SNPs 
      // 3) Haplotypes
      // 4) Covariates
      // 5) Interactions of the above
      // 6) QFAM variables

      // Populate this vector with terms for this
      // individual
           
      skip = false;
      
      vector_t trow(np);
      
      for (int p = 0; p < np; p++)
	{
	  
	  int pType = type[p];
	  
	  switch ( pType )
	    {
	    case INTERCEPT :  
	      trow[p] = buildIntercept();
	      break;
	    case ADDITIVE :
	      trow[p] = buildAdditive( person, order[p] );
	      break;
	    case DOMDEV :
	      trow[p] = buildDominance( person, order[p] );
	      break;
	    case HAPLOTYPE :
	      trow[p] = buildHaplotype(i, order[p] );
	      break;
	    case SEX :
	      trow[p] = buildSex(person);
	      break;
	    case COVARIATE :
	      trow[p] = buildCovariate( person, order[p] );
	      break;
	    case INTERACTION :
	      trow[p] = buildInteraction( person, order[p], trow );
	      break;
	    case QFAM :
	      trow[p] = buildQFAM( person );
	      break;

	    }
	}

      if (skip) 
	{
	  miss[i] = true;
	  skip = false;
	  continue;
	}
            
	
      ////////////////////////////
      // Add row to design matrix
      
      X.push_back(trow);
            
    }
      

  /////////////////////////////////////////////////
  // Set number of non-missing individuals

  nind = X.size();


  //////////////////////////////////
  // Apply a parameter list filter?
  
  if ( par::glm_user_parameters )
    {

      // Intercept always fixed in
      
      np = 1;
      vector<string> label2 = label;
      label.clear();
      label.push_back("M");
      int np2 = label2.size();
      
      for (int i=0; i< par::parameter_list.size(); i++)
	{
	  if ( par::parameter_list[i] >= 1 
	       && par::parameter_list[i] < np2 )
	    {
	      np++;
	      label.push_back(label2[ par::parameter_list[i] ]);
	    }
	}

     
      // For each individual

      for ( int i = 0 ; i < X.size() ; i++)
	{

	  vector_t X2(1);
	  
	  X2[0] = 1;
	  
	  for ( int j = 0 ; j < par::parameter_list.size() ; j++)
	    {
	    if ( par::parameter_list[j] >= 1 
		 && par::parameter_list[j] < np2 )
	      {
		X2.push_back(X[i][ par::parameter_list[j] ]);
	      }
	    }
	  X[i] = X2;
	}
    }


  /////////////////////////////////////////
  // VIF-based check for multicollinearity

  all_valid = checkVIF();


  ///////////////////////
  // Add Y variable also

  setDependent();

  // Now we are ready to perform the analysis

  if (par::verbose)
    {
      cout << "X design matrix\n";
      display(X);
      cout << "\n";
    }
      
}


vector<bool> Model::validParameters()
{  

  // Empty model?
  if (np==0 || nind==0) 
    {
      vector<bool> v(np,false);
      all_valid = false;
      return v;
    }
  

  // Display covariance matrix in verbose mode

  if ( par::verbose )
    {
      cout << "Covariance matrix of estimates\n";
      display(S);
      cout << "\n";
    }

  // Check for multicollinearity
  
  // For each term, see that estimate is not too strongly (r>0.99)
  // correlated with another, starting at last

  valid.resize(np);
  
  for (int i = 1; i<np; i++)
    {
      valid[i] = true;
      if ( S[i][i] < 1e-20 ) { valid[i] = all_valid = false; } 
      else if ( ! realnum(S[i][i]) ) { valid[i] = all_valid = false; }
    }

  if ( all_valid ) 
    for (int i = np-1; i>0; i--)
      {
	for (int j = i-1; j>=0; j--)
	  {	  
	    if ( S[i][j] / sqrt( S[i][i] * S[j][j] ) > 0.99999 )
	      {
		valid[i] = false;
		all_valid = false;
		break;
	      }
	  }
      }
  return valid;
}


double Model::getStatistic()
{
  if (all_valid)
    {
      return ( coef[testParameter] * coef[testParameter] ) 
	/ S[testParameter][testParameter];
    }
  else return 0;
}


// **********************************************
// *** Function to aid testing linear 
// *** hypotheses after estimating of a
// *** regression
// **********************************************

double Model::linearHypothesis(matrix_t & H, vector_t & h)
{
  
  // (H v - h)' (H S H')^-1 (H b - h) ~ X^2 with j df 
  // where H = constraint matrix (j x (p+1)
  //       h = null coefficient values
  //       S = estimated covariance matrix of coefficients

  //   return ( H*b - h ).transpose() 
  //     * ( H*v*H.transpose() ).inverse()
  //     * ( H*b - h ) ;
  
  int nc = h.size(); // # of constraints
  
  // 1. Calculate Hb-h  

  vector_t outer;
  outer.resize(nc,0);
  
  for (int r = 0; r < nc; r++)
    for (int c = 0; c < np; c++)
      outer[r] += H[r][c] * coef[c];
  
  for (int r = 0; r < nc; r++)
    outer[r] -= h[r];
  
  // 2. Calculate HVH'
  
  matrix_t tmp;
  sizeMatrix(tmp,nc,np);
  
  for (int r = 0; r < nc; r++)
    for (int c = 0; c < np; c++)
      for (int k = 0; k < np; k++)
	tmp[r][c] += H[r][k] * S[k][c];
  
  matrix_t inner;
  sizeMatrix(inner,nc,nc);
  
  for (int r = 0; r < nc; r++)
    for (int c = 0; c < nc; c++)
      for (int k = 0; k < np; k++)
	inner[r][c] += tmp[r][k] * H[c][k];

  bool flag = true;
  inner = svd_inverse(inner,flag);
  if ( !flag ) all_valid = false;

  vector_t tmp2;
  tmp2.resize(nc,0);
  
  for (int c = 0; c < nc; c++)
    for (int k = 0; k < nc; k++)
      tmp2[c] += outer[k] * inner[k][c];

  double result = 0;
  
  for (int r = 0; r < nc; r++)
    result += tmp2[r] * outer[r];

  return result;
    
}


bool Model::checkVIF()
{

  // Calculate correlation matrix for X
  // Skip intercept

  int p = X.size();
  if (p<2) return false;

  int q = X[0].size() - 1; 	
  if ( q < 2 ) return true;

  vector_t m(q);
  matrix_t c;
  sizeMatrix(c,q,q);

  for (int i=0; i<p; i++)
    for (int j=0; j<q; j++)
      m[j] += X[i][j+1];

  for (int j=0; j<q; j++)
    m[j] /= (double)p;
  
  for (int i=0; i<p; i++)
    for (int j1=0; j1<q; j1++)
      for (int j2=j1; j2<q; j2++)
	c[j1][j2] += ( X[i][j1+1] - m[j1] ) * ( X[i][j2+1] - m[j2] );  
  
  for (int j1=0; j1<q; j1++)
    for (int j2=j1; j2<q; j2++)
      c[j1][j2] /= (double)(p-1);
  
  for (int j1=0; j1<q; j1++)
    for (int j2=j1+1; j2<q; j2++)
      {
	c[j1][j2] /= sqrt( c[j1][j1] * c[j2][j2] );
	c[j2][j1] = c[j1][j2];
	
 	if ( c[j2][j1] > 0.999 ) 
 	  {
	    if (par::verbose)
	      cout << "individual element > 0.999, skipping VIF\n";
	    return false;
	  }
      }		
  

  // Any item with zero variance?

  for (int j=0; j<q; j++)
    {
      if ( c[j][j] == 0 || ! realnum( c[j][j] ) )
	return false;
      c[j][j] = 1;
    }
  
  // Get inverse
  bool flag = true;
  c = svd_inverse(c,flag);
  if ( ! flag ) all_valid = false;

  if (par::verbose)
    {
      cout << "VIF on diagonals\n";
      display(c);
      cout << "\n";
    }

  // Calculate VIFs
  double maxVIF = 0;
  for (int j=0;j<q;j++)
    {
      
      // r^2 = 1 - 1/x where x is diagonal element of inverted
      // correlation matrix
      // As VIF = 1 / ( 1 - r^2 ) , implies VIF = x

      if ( c[j][j] > par::vif_threshold )
	return false;
    }
  
  
  return true;
  
}




double Model::buildIntercept()
{
  return 1;
}

double Model::buildAdditive(Individual * person , int snp )
{

  // Additive effects (assuming individual-major mode)
      
  int s = additive[snp];  

  bool i1 = person->one[s];
  bool i2 = person->two[s];
      
  if ( xchr[snp] )
    {

      /////////////////////////
      // X chromosome coding
      
      if ( person->sex ) // male
	{
	  if ( i1 ) 
	    {
	      if ( ! i2 ) 
		{
		  skip = true;
		  return 0;
		}
	      else
		return mA; 	      
	    }
	  else
	    {
	      if ( i2 )
		{
		  // This should not happen...
		  skip = true;
		  return 0;
		}		  
	      else
		return mB; 
	    }
	    }
      else // female x-chromosome
	{
	  if ( i1 ) 
	    {
	      if ( ! i2 ) 
		{
		  skip = true;
		  return 0;
		}
	      else
		return mAA; 	      
	    }
	  else
	    {
	      if ( i2 )
		return mAB; // het 
	      else
		return mBB; // hom
	    }		  
	}
      
    }
  else if ( haploid[snp] )
    {
      
      ///////////////////
      // Haploid coding
      
      if ( i1 ) 
	{
	  if ( ! i2 ) 
	    {
	      skip = true;
	      return 0;
	    }
	  else
	    return 0; 	      
	}
      else
	{
	  if ( i2 )
	    {
	      // haploid het
	      skip = true;
	      return 0;
	    }		  
	  else
	    return 1;
	    }
      
    }
  else 
    {
      
      ///////////////////////
      // Autosomal coding
      
      if ( i1 ) 
	{
	  if ( ! i2 ) 
	    {
	      skip = true;
	      return 0;
	    }
	  else
	    return mAA;
	}
      else
	{
	  if ( i2 )
	    return mAB; // het 
	  else
	    return mBB; // hom
	}
      
    }
  
}


double Model::buildDominance(Individual * person, int snp)
{

  ////////////////////
  // Dominance effects
  
  int s = dominance[snp];
  bool i1 = person->one[s];
  bool i2 = person->two[s];
	  
  if ( i1 ) 
    {
      if ( ! i2 ) 
	{
	  skip = true;
	  return 0;
	}
      else
	return 0;
    }
  else
    {
      if ( i2 )
	return 1; // het 
      else
	return 0; // hom
    }
  
}

double Model::buildHaplotype(int i, int h )
{
  
  ////////////////////
  // Haplotype dosage

  if ( P->haplo->include[i] )
    return P->haplo->dosage(i,haplotype[h]);

  // No valid haplotypes for this person
  skip = true;
  return 0;
    
}

double Model::buildSex(Individual * person )
{
      
  ////////////////////////////////////
  // Sex effect (automatically set for 
  // X chromosome models)
  
  if ( person->sex )
    return 1;
  else 
    return 0;
  
}

double Model::buildCovariate(Individual * person, int j)
{

  /////////////
  // Covariates
  
  return person->clist[covariate[j]];
      
}

double Model::buildInteraction(Individual * person, int j, vector_t & trow )
{
  ///////////////
  // Interactions
  
  return trow[ interaction[j].p1 ] * trow[ interaction[j].p2 ];
    
}

double Model::buildQFAM(Individual * person)
{

  ///////////////
  // QFAM
  
  if ( par::QFAM_total )
    return person->T;
  else if ( par::QFAM_between )
    return person->family->B;
  else if ( par::QFAM_within1 || par::QFAM_within2 ) 
    return person->W;
  else
    error("Internal problem with QFAM model specification");
}



void Model::setCluster()
{
  cluster=true;
  clst.clear();
  nc=0;

  map<int,int> novelc;
  
  for (int i=0; i<P->n; i++)
    if (!miss[i])
      {
	int c = P->sample[i]->sol;
	map<int,int>::iterator m = novelc.find(c);
	if ( m == novelc.end() )
	  {
	    clst.push_back( nc );
	    novelc.insert( make_pair( c, nc ) );
	    ++nc;
	  }
	else
	  {
	    clst.push_back( m->second );
	  }
      }
  
  // We need at least two clusters:
  if ( novelc.size() == 1 )
    noCluster();
}

void Model::noCluster()
{
  cluster=false;
  clst.clear();
  nc=0;
}

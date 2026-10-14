

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#include <iostream>
#include <iomanip>
#include <cmath>

#include "logistic.h"
#include "plink.h"
#include "helper.h"
#include "options.h"
#include "stats.h"

LogisticModel::LogisticModel(Plink * p_)
{
  P = p_;
  nc = 0;
  cluster = false;
}

void LogisticModel::setDependent()
{

  // Set phenotype to 'aff' variable
  Y.clear();

  for (int i=0; i<P->n; i++)
    {
      if ( !miss[i] )
	{
	  if ( P->sample[i]->pperson->aff ) 
	    Y.push_back( 1 ) ;
	  else
	    Y.push_back( 0 ) ;
	}
    }
  
  nind = Y.size();
  
  p.resize(nind);
  V.resize(nind);
  
}


void LogisticModel::pruneY()
{

  //////////////////////////////////
  // Prune out rows that are missing

  if ( miss.size() != Y.size() ) 
    error("Internal error: bad call to Model::pruneY()");

  vector<int> Y2;
  for (int i=0; i<Y.size(); i++)
    if ( ! miss[i] ) 
      Y2.push_back(Y[i]);
  Y = Y2;
}

void LogisticModel::fitLM() 
{

  coef.resize(np);
  sizeMatrix(S,np,np);

  if (np==0 || nind==0 || ! all_valid )
    return;

  if (par::verbose)
    {
      for (int i=0; i<nind; i++)
	{
	  cout << i << "\t"
	       << Y[i] << "\t";
	  for (int j=0; j<np; j++)
	    cout << X[i][j] << "\t";
	  cout << "\n";
	}
    }

  ///////////////////////////////////////
  // Newton-Raphson to fit logistic model
    
    bool converge = false;
    int it = 0;

    while ( ! converge && it < 20 ) 
      {
	
	// Determine p and V
	for (int i=0; i<nind; i++)
	  {
	    double t = 0;
	    for (int j=0; j<np; j++)
	      t += coef[j] * X[i][j];	    
	    p[i] = 1/(1+exp(-t));
	    V[i] = p[i] * (1-p[i]);
	  }
	
	// Update coefficients
	// b <- b +  solve( t(X) %*% V %*% X ) %*% t(X) %*% ( y - p ) 

	matrix_t T;
	sizeMatrix(T,np,np);

	for (int j=0; j<np; j++)
	  for (int k=j; k<np; k++) 
	    {
	      double sum = 0;
	      for (int i=0; i<nind; i++)
		sum += X[i][j] * V[i] * X[i][k] ;
	      T[j][k] = T[k][j] = sum;	      
	    }

	bool flag = true;
	T = svd_inverse(T,flag);
	if ( ! flag ) 
	  {
	    all_valid = false;
	    return;
	  }

	matrix_t T2;
	// Resize and set elements to 0
	sizeMatrix(T2,np,nind); 
	
	// note implicit transpose of X
	for (int i=0; i<np; i++)
	  for (int j=0; j<nind; j++)
	    for (int k=0; k<np; k++)
	      T2[i][j] += T[i][k] * X[j][k];  
		
	vector_t t3(nind);
	for (int i=0; i<nind; i++) 
	  t3[i] = Y[i] - p[i];
	
	vector_t ncoef(np);
	for (int j=0; j<np; j++) 
	  for (int i=0; i<nind; i++) 
	    ncoef[j] += T2[j][i] * t3[i];

	// Update coefficients, and check for 
	// convergence
	double delta = 0;
	for (int j=0; j<np; j++) 	
	  {
	    delta += abs(ncoef[j]);
	    coef[j] += ncoef[j];
	  }

	if ( delta < 1e-6 )
	  converge = true;

	// Next iteration
	it++;
      }


    /////////////////////////////////////////
    // Obtain covariance matrix of estimates
      
    // S <- solve( t(X) %*% V %*% X )    
      
    // Transpose X and multiple by diagonal V
    matrix_t Xt;
    sizeMatrix(Xt, np, nind);
    for (int i=0; i<nind; i++)
      for (int j=0; j<np; j++) 
	Xt[j][i] = X[i][j] * V[i];
    
    multMatrix(Xt,X,S);  
    bool flag = true;
    S = svd_inverse(S,flag);     
    if ( ! flag ) 
      {
	all_valid = false;
	return;
      }
    if ( cluster ) 
      HuberWhite();

    if (par::verbose)
      {
	cout << "beta\n";
	display(coef);
	cout << "Sigma\n";
	display(S);
	cout << "\n";
      }
}

vector_t LogisticModel::getCoefs()
{
  return coef;
}

vector_t LogisticModel::getVar()
{  
  vector_t var(np);
  for (int i=0; i<np; i++) 
    var[i] = S[i][i];
  return var;
}

vector_t LogisticModel::getSE()
{  
  vector_t var(np);
  for (int i=0; i<np; i++) 
    var[i] = sqrt(S[i][i]);
  return var;
}

void LogisticModel::reset()
{
  np=0;
  nind=0;
  coef.clear();
  S.clear();
  Y.clear();
  X.clear();
  miss.clear();
}

void LogisticModel::displayResults(ofstream & OUT, Locus * loc)
{

  vector_t var;
  
  if ( all_valid ) 
    var = getVar();
  else
    {
      var.clear();
      var.resize(np,0);
    }
  

  for (int p=1; p<np; p++) // Skip intercept
    {
 
      bool okay = var[p] < 1e-20 || !realnum(var[p]) ? false : all_valid;
      
      double se = 0; 
      double Z = 0;
      double pvalue = 1;
      
      if (okay)
	{
	  se = sqrt(var[p]);
	  Z = coef[p] / se;
	  //	  pvalue = pT(Z,Y.size()-np);
	  pvalue = chiprobP(Z*Z,1);
	}
      
      // If filtering p-values
      if ( (!par::pfilter) || pvalue <= par::pfvalue ) 
	{	 

	  // Skip covariates?
	  if ( par::no_show_covar && p != testParameter )
	    continue;

	  OUT << setw(4) << loc->chr  << " " 
	      << setw(par::pp_maxsnp) << loc->name << " " 
	      << setw(10) << loc->bp << " " 
	      << setw(4) << loc->allele1 << " "
	      << setw(10) << label[p] << " "
	      << setw(8) << Y.size() << " ";
	  
	  if (okay)
	    {
	      if ( par::return_beta )
		OUT << setw(10) << coef[p] << " ";
	      else
		OUT << setw(10) << exp(coef[p]) << " ";
	      
	      if (par::display_ci)
		{
		  OUT << setw(8) << se << " ";

		  if ( par::return_beta )
		    OUT << setw(8) << coef[p] - par::ci_zt * se << " "
			<< setw(8) << coef[p] + par::ci_zt * se << " ";	    
		  else
		    OUT << setw(8) << exp(coef[p] - par::ci_zt * se) << " "
			<< setw(8) << exp(coef[p] + par::ci_zt * se) << " ";	    

		}	      

	      OUT << setw(12) << Z << " "
		  << setw(12) << pvalue;
	    }
	  else
	    {
	      OUT << setw(10) << "NA" << " ";
	      
	      if (par::display_ci)
		OUT << setw(8) << "NA" << " "  
		    << setw(8) << "NA" << " "
		    << setw(8) << "NA" << " ";	    
	      
	      OUT << setw(12) << "NA" << " "
		  << setw(12) << "NA";
	    }
	  
	  OUT << "\n";	
	}
    }
}


double LogisticModel::getPValue()
{  
  vector_t var = getVar();
  bool okay = var[testParameter] < 1e-20 || !realnum(var[testParameter]) ? false : all_valid;

  if (all_valid)
    {
      double se = sqrt(var[testParameter]);
      double Z = coef[testParameter] / se;	  
      return chiprobP(Z*Z,1);
    }
  else return 1;
}

vector_t LogisticModel::getPVals()
{
  int tmp = testParameter;
  vector_t res;
  for ( testParameter = 1; testParameter < np; testParameter++)
    res.push_back( getPValue() );
  testParameter = tmp;
  return res;
}


double LogisticModel::getLnLk()
{
  // Return -2 * sample log-likelihood
  // We assume the model is fit, and all Y's are either 0 or 1

  double lnlk = 0;
  
  for (int i=0; i<nind; i++)
    {
      double t = 0;
      for (int j=0; j<np; j++)
	t += coef[j] * X[i][j];	    
      lnlk += Y[i] == 1 ? log( 1/(1+exp(-t))) : log(1 - (1/(1+exp(-t))) );
    }
  
  return -2 * lnlk;
  
}


void LogisticModel::HuberWhite()
{
  
  // Calculate sandwich variance estimators, potentially allowing for
  // clustered data
  
  // Works to update the S matrix, variance/covariance matrix
  
  // Originally, S will contain this, uncorrected

  // Calcuate S = (XtX)^-1
  
  
  matrix_t S0 = S; 

  vector<vector_t> sc(nc);
  for (int i=0; i<nc; i++)
    sc[i].resize(np,0);

  for (int i=0; i<nind; i++)
    {
      double err = Y[i] - p[i];      
      for (int j=0; j<np; j++)
	sc[clst[i]][j] += err * X[i][j];
    }
  
  matrix_t meat;
  sizeMatrix(meat, np, np);
  for (int k=0; k<nc; k++)
    {      
      for (int i=0; i<np; i++)
 	for (int j=0; j<np; j++)
 	  meat[i][j] += sc[k][i] * sc[k][j];
      
    }

  matrix_t tmp1;
  multMatrix( S0 , meat, tmp1);
  multMatrix( tmp1 , S0, S);
  
}

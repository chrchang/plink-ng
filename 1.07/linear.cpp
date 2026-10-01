

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

#include "linear.h"
#include "helper.h"
#include "options.h"
#include "stats.h"

LinearModel::LinearModel(Plink * p_)
{
  P = p_;

  nc = 0;
  cluster = false;
  RSS = -1;
}

void LinearModel::setDependent()
{
  // Set dependent variable and intercept
  Y.clear();

  for (int i=0; i<P->n; i++)
    if (!miss[i])
      Y.push_back( P->sample[i]->pperson->phenotype ) ;

}

void LinearModel::pruneY()
{

  //////////////////////////////////
  // Prune out rows that are missing
  
  vector<double> Y2;
  for (int i=0; i<Y.size(); i++)
    if ( ! miss[i] ) 
      Y2.push_back(Y[i]);
  Y = Y2;
}

void covsrt(matrix_t & covar, vector<bool> &ia, const int mfit)
{
  int i,j,k;
  
  int ma=ia.size();
  for (i=mfit;i<ma;i++)
    for (j=0;j<i+1;j++) covar[i][j]=covar[j][i]=0.0;
  k=mfit-1;
  for (j=ma-1;j>=0;j--) {
    if (ia[j]) {
      for (i=0;i<ma;i++) SWAP(covar[i][k],covar[i][j]);
      for (i=0;i<ma;i++) SWAP(covar[k][i],covar[j][i]);
      k--;
    }
  }
}

void gaussj(matrix_t & a, matrix_t & b)
{
  int i,icol,irow,j,k,l,ll;
  double big,dum,pivinv;
  
  int n=a.size();
  int m=b[0].size();
  vector_t indxc(n),indxr(n),ipiv(n);
  for (j=0;j<n;j++) ipiv[j]=0;
  for (i=0;i<n;i++) {
    big=0.0;
    for (j=0;j<n;j++)
      if (ipiv[j] != 1)
	for (k=0;k<n;k++) {
	  if (ipiv[k] == 0) {
	    if (fabs(a[j][k]) >= big) {
	      big=fabs(a[j][k]);
	      irow=j;
	      icol=k;
	    }
	  }
	}
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l]);
      for (l=0;l<m;l++) SWAP(b[irow][l],b[icol][l]);
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol][icol] == 0.0) error("gaussj: Singular Matrix");
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=0;l<n;l++) a[icol][l] *= pivinv;
    for (l=0;l<m;l++) b[icol][l] *= pivinv;
    for (ll=0;ll<n;ll++)
      if (ll != icol) {
	dum=a[ll][icol];
	a[ll][icol]=0.0;
	for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
	for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
      }
  }
  for (l=n-1;l>=0;l--) {
    if (indxr[l] != indxc[l])
      for (k=0;k<n;k++)
	SWAP(a[k][(int)indxr[l]],a[k][(int)indxc[l]]);
  }
}

void lfit(vector_t &x, vector_t &y, vector_t &sig, vector_t &a,
 	  vector<bool> &ia, matrix_t &covar, double &chisq, matrix_t & X)
{
  int i,j,k,l,m,mfit=0;
  double ym,wt,sum,sig2i;
  
  int ndat=x.size();
  int ma=a.size();
  vector_t afunc(ma);
  matrix_t beta;
  sizeMatrix(beta,ma,1);
  for (j=0;j<ma;j++)
    if (ia[j]) mfit++;
  if (mfit == 0) error("lfit: no parameters to be fitted");
  for (j=0;j<mfit;j++) {
    for (k=0;k<mfit;k++) covar[j][k]=0.0;
    beta[j][0]=0.0;
  }
  for (i=0;i<ndat;i++) {
    afunc = X[i];

    ym=y[i];
    if (mfit < ma) {
      for (j=0;j<ma;j++)
	if (!ia[j]) ym -= a[j]*afunc[j];
    }
    sig2i=1.0/SQR(sig[i]);
    for (j=0,l=0;l<ma;l++) {
      if (ia[l]) {
	wt=afunc[l]*sig2i;
	for (k=0,m=0;m<=l;m++)
	  if (ia[m]) covar[j][k++] += wt*afunc[m];
	beta[j++][0] += ym*wt;
      }
    }
  }
  for (j=1;j<mfit;j++)
    for (k=0;k<j;k++)
      covar[k][j]=covar[j][k];
  vector<vector<double> >  temp;
  sizeMatrix(temp,mfit,mfit);
  for (j=0;j<mfit;j++)
    for (k=0;k<mfit;k++)
      temp[j][k]=covar[j][k];
  gaussj(temp,beta);
  for (j=0;j<mfit;j++)
    for (k=0;k<mfit;k++)
      covar[j][k]=temp[j][k];
  for (j=0,l=0;l<ma;l++)
    if (ia[l]) a[l]=beta[j++][0];
  chisq=0.0;
  for (i=0;i<ndat;i++) {
    afunc = X[i];
    sum=0.0;
    for (j=0;j<ma;j++) sum += a[j]*afunc[j];
    chisq += SQR((y[i]-sum)/sig[i]);
  }
  covsrt(covar,ia,mfit);
}


void LinearModel::setVariance()
{
  varY=0;
  meanY = 0;
  int actualN=0;

  for (int i=0; i<nind; i++)
    {
      actualN++;
      meanY += Y[i];
    }
  if (actualN==0)
    {
      varY=0;
      return;
    }
  meanY /= (double)actualN;
  
  for (int i=0; i<nind; i++)
    {
      varY += (Y[i] - meanY) * (Y[i] - meanY);	
    }
  varY /= (double)(actualN-1);

  if ( actualN != nind ) 
    error("actualN <> nind...");
}

void LinearModel::standardise()
{
  
  // Get mean and variance for all variable
  double sdY = sqrt(varY);
  for (int i=0; i<nind; i++)
    Y[i] = ( Y[i] - meanY ) / sdY;

  vector_t mX(np,0);
  vector_t sdX(np,0);

  // Standardise all predictors, except intercept
  for (int i=0; i<nind; i++)
    for (int j=1; j<np; j++)
      mX[j] += X[i][j];
  for (int j=1; j<np; j++)
    mX[j] /= nind;
  for (int i=0; i<nind; i++)
    for (int j=1; j<np; j++)
      sdX[j] += (X[i][j] - mX[j]) * (X[i][j] - mX[j]);
  for (int j=1; j<np; j++)
    {
      sdX[j] = sqrt(sdX[j]/(nind-1));
      if ( sdX[j]==0 )
	sdX[j] = 1;
    }
  for (int i=0; i<nind; i++)
    for (int j=1; j<np; j++)
      X[i][j] = ( X[i][j] - mX[j] ) / sdX[j];
}


void LinearModel::fitLM() 
{

  if (par::verbose)
    {
      for (int i=0; i<nind; i++)
	{
	  cout << "VO " << i << "\t"
	       << Y[i] << "\t";
 	  for (int j=0; j<np; j++)
	    cout << X[i][j] << "\t";
	  cout << "\n";
	}
    }

//   cout << "LM VIEW\n";
//       display(Y);
//       display(X);
//       cout << "---\n";

  coef.resize(np);
  sizeMatrix(S,np,np);

  if ( np==0 || nind==0 || ! all_valid )
    {
      return;
    }
  
  setVariance();

  if ( par::standard_beta )
    standardise();


  sig.resize(nind, sqrt(1.0/sqrt((double)nind)) );
  
  w.resize(np);
  sizeMatrix(u,nind,np);
  sizeMatrix(v,np,np);
  

  //  Perform "svdfit(C,Y,sig,b,u,v,w,chisq,function)"
  
  int i,j;
  const double TOL=1.0e-13;
  double wmax,tmp,thresh,sum;
    
  vector_t b(nind),afunc(np);
  for (i=0;i<nind;i++) {
    afunc = X[i];
    tmp=1.0/sig[i];
    for (j=0;j<np;j++) 
      u[i][j]=afunc[j]*tmp;     
    b[i]=Y[i]*tmp;
  }

  bool flag = svdcmp(u,w,v);
  
  if ( ! flag ) 
    {
      all_valid = false;
      return;
    }

  wmax=0.0;
  for (j=0;j<np;j++)
    if (w[j] > wmax) wmax=w[j];
  thresh=TOL*wmax;
  for (j=0;j<np;j++)
    if (w[j] < thresh) w[j]=0.0;
  
  svbksb(u,w,v,b,coef);
  
  chisq=0.0;
  for (i=0;i<nind;i++) {
    afunc=X[i];
    sum=0.0;
    for (j=0;j<np;j++) sum += coef[j]*afunc[j];
    chisq += (tmp=(Y[i]-sum)/sig[i],tmp*tmp);
  }



  /////////////////////////////////////////
  // Obtain covariance matrix of estimates


  // Robust cluster variance estimator
  // V_cluster = (X'X)^-1 * \sum_{j=1}^{n_C} u_{j}' * u_j * (X'X)^-1 
  // where u_j = \sum_j cluster e_i * x_i 

  // Above, e_i is the residual for the ith observation and x_i is a
  // row vector of predictors including the constant.

  // For simplicity, I omitted the multipliers (which are close to 1)
  // from the formulas for Vrob and Vclusters.

  // The formula for the clustered estimator is simply that of the
  // robust (unclustered) estimator with the individual ei*s replaced
  // by their sums over each cluster. 

  // http://www.stata.com/support/faqs/stat/cluster.html
  // SEE http://aje.oxfordjournals.org/cgi/content/full/kwm223v1#APP1

  // Williams, R. L. 2000.  A note on robust variance estimation for
  // cluster-correlated data. Biometrics 56: 64

  //  t ( y - yhat X  ) %*%  ( y - yhat)  / nind - np
  // = variance of residuals 
  // j <- ( t( y- m %*% t(b) ) %*% ( y - m %*% t(b) ) ) / ( N - p ) 
  // print( sqrt(kronecker( solve( t(m) %*% m ) , j )  ))
  

  ////////////////////////////////////////////////
  // OLS variance estimator = s^2 * ( X'X )^-1
  // where s^2 = (1/(N-k)) \sum_i=1^N e_i^2
  
  // 1. Calcuate S = (X'X)^-1
  
  matrix_t Xt;
  sizeMatrix(Xt, np, nind);
  for (int i=0; i<nind; i++)
    for (int j=0; j<np; j++) 
      Xt[j][i] = X[i][j];
  matrix_t S0; 
  multMatrix(Xt,X,S0);
  flag = true;
  S0 = svd_inverse(S0,flag);  
  if ( ! flag ) 
    {
      all_valid = false;
      return;
    }

  if (par::verbose)
    {
      cout << "beta...\n";
      display(coef);
      cout << "Sigma(S0b)\n";
      display(S0);
      cout << "\n";
    }


  ////////////////////////
  // Calculate s^2 (sigma)

  if (!cluster)
    {
      double sigma= 0.0;
      for (int i=0; i<nind; i++)
	{
	  double partial = 0.0;
	  for (int j=0; j<np; j++)
	    partial += coef[j] * X[i][j];
	  partial -= Y[i];
	  sigma += partial * partial;
	}
      sigma /= nind-np;	
            
      for (int i=0; i<np; i++)
	for (int j=0; j<np; j++)
	  S[i][j] = S0[i][j] * sigma;
    }
  

  ///////////////////////////
  // Robust-cluster variance

  if (cluster)
  {
    
    vector<vector_t> sc(nc);
    for (int i=0; i<nc; i++)
      sc[i].resize(np,0);

    for (int i=0; i<nind; i++)
      {
	double partial = 0.0;
	for (int j=0; j<np; j++)
	  partial += coef[j] * X[i][j];
	partial -= Y[i];
	
	for (int j=0; j<np; j++)
	  sc[clst[i]][j] += partial * X[i][j];
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

  
  if (par::verbose)
    {
      cout << "coefficients:\n";
      display(coef);
      cout << "var-cov matrix:\n";
      display(S);
      cout << "\n";
    }
  
}


void LinearModel::fitUnivariateLM() 
{

  if (par::verbose)
    {
      cout << "LM VIEW\n";
      display(Y);
      display(X);
      cout << "---\n";
    }

  // Speed-up version for univariate case Has set set coef and S
  
  if (np!=2 || nind==0)
    return;
  
  coef.resize(2);
  sizeMatrix(S,2,2);
  
  double x_mean=0, x_var=0;
  double y_mean=0, y_var=0;
  double y_x_covar=0;
  

  /////////////////////////////
  // Iterate over individuals

  // X and Y
  
  for (int i=0; i<nind; i++)
    {
      y_mean += Y[i];
      x_mean += X[i][1];
    }

  x_mean /= (double)nind;
  y_mean /= (double)nind;      
  
      
  for (int i=0; i<nind; i++)
    {
      double ty = Y[i] - y_mean;
      double tx = X[i][1] - x_mean;
      y_var += ty*ty;
      x_var += tx*tx;
      y_x_covar += tx * ty;
    }

      
  y_var /= (double)nind - 1;
  x_var /= (double)nind - 1;
  y_x_covar /= (double)nind - 1;
  
  
  // Do not set intercept; only the univariate coefficient
  coef[1] = y_x_covar / x_var;
  S[1][1] = (y_var/x_var - (y_x_covar*y_x_covar)/(x_var*x_var) ) / (nind-2);
  
}


vector_t LinearModel::getCoefs()
{
 return coef;
}


vector_t LinearModel::getVar()
{

  double multiplier = 1;

  if (cluster) 
    multiplier = (double)(nc)/((double)(nc-np));

  multiplier = 1;

  //cout << "Mult = " << multiplier << "\n";  

  vector_t var(np);
  for (int i=0; i<np; i++) 
   var[i] = multiplier * S[i][i];

  return var;
}


vector_t LinearModel::getSE()
{
  
  double multiplier = 1;
  
  if (cluster) 
    multiplier = (double)(nc)/((double)(nc-np));

  multiplier = 1;

  vector_t var(np);
  for (int i=0; i<np; i++) 
    var[i] = sqrt ( multiplier * S[i][i] ) ;
  
  return var;
}


void LinearModel::reset()
{
  np=0;
  nind=0;
  coef.clear();
  w.clear();
  S.clear();
  Y.clear();
  X.clear();
  miss.clear();
}


void LinearModel::displayResults(ofstream & OUT, Locus * loc)
{
  
  vector_t var;
  
  if ( all_valid ) 
    var = getVar();
  else
    {
      var.clear();
      var.resize(np,0);
    }
  
  for (int p=1; p<np; p++) // skip intercept
    {
      
      bool okay = var[p] < 1e-20 || !realnum(var[p]) ? false : all_valid;
      
      double se = 0; 
      double Z = 0;
      double pvalue = 1;
      
      if (okay)
	{
	  se = sqrt(var[p]);
	  Z = coef[p] / se;
	  pvalue = pT(Z,Y.size()-np);
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
	      OUT << setw(10) << coef[p] << " ";
	      
	      if (par::display_ci)
		
		OUT << setw(8) << se << " "
		    << setw(8) << coef[p] - par::ci_zt * se << " "
		    << setw(8) << coef[p] + par::ci_zt * se << " ";	    
	      
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


double LinearModel::calculateRSS()
{

  // Might already be calculated?

  if ( RSS >= 0 ) return RSS;

  // Calculate residual sum of squares (RSS)
  
  RSS = 0;
  
  for (int i=0; i<nind; i++)
    {

      double t = Y[i];
      
      for ( int p=0; p<np; p++)
	t -= coef[p] * X[i][p];
      
      t *= t;
      RSS += t; 
    }
  
  return RSS;
}


double LinearModel::calculateRSquared()
{

  // Return coefficient of determination. Ifnot already calculated,
  // get residual sum of squares first (set to -1)
  
  if ( RSS < 0 ) 
    RSS = calculateRSS();
  
  double SSy = varY * (nind-1);

  double r = ( SSy - RSS ) / SSy;

  return r > 0 ? ( r > 1 ? 1 : r ) : 0;

}

double LinearModel::calculateAdjustedRSquared()
{
  
  double ra =  1 - ( (double)(nind-1)/(double)(nind-np-1) ) 
    * ( 1 - calculateRSquared() );

  return ra > 0 ? ( ra > 1 ? 1 : ra ) : 0; 

}

double LinearModel::calculateMallowC(LinearModel * submodel)
{

  // Mallow's C = RSSm / S^2 + 2(m+1)-n
  // where S^2 = RSSk / (n-k-1);
  
  double Sk = calculateRSS() / ( nind - np - 1); 
  return ( submodel->calculateRSS() / Sk ) + 2 * ( submodel->np+1)-nind;
}


double LinearModel::calculateFTest(LinearModel * submodel)
{

  double RSSk = calculateRSS();
  double RSSm = submodel->calculateRSS();
  
  return ( ( RSSm - RSSk ) / (double)( np - submodel->np ) )
    / ( RSSk / (double)(nind - np - 1 ) );
}

double LinearModel::getPValue()
{  
  vector_t var = getVar();
  bool okay = var[testParameter] < 1e-20 || !realnum(var[testParameter]) ? false : all_valid;

  if (all_valid)
    {
      double se = sqrt(var[testParameter]);
      double Z = coef[testParameter] / se;	  
      return pT(Z,Y.size()-np);
    }
  else return 1;
}

vector_t LinearModel::getPVals()
{
  int tmp = testParameter;
  vector_t res;
  for (testParameter = 1; testParameter < np; testParameter++)
    res.push_back( getPValue() );
  testParameter = tmp;
  return res;
}

void LinearModel::HuberWhite()
{
  //
}

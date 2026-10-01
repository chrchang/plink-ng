



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
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <functional>
#include <cmath>

#include "plink.h"
#include "options.h"
#include "sets.h"
#include "helper.h"
#include "stats.h"

#include "crandom.h"
#include "linear.h"
#include "logistic.h"

typedef vector<long double> vector_tld;

// Function that implements Pillai's (1964) approximation to upper
// distribution of the largest canonical correlatio

// Function (bico) to estimate combinations choose(n,k)
double factln(int n)
{
    double gammln(double xx);
    static double a[101];

    if (n < 0) error("Negative factorial in routine factln");
    if (n <= 1) return 0.0;
    if (n <= 100 ) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
    else return gammln(n+1.0);
}
double bico(int n, int k)
{
    double factln(int n);
    return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}

// Beta function
double betaln(double z, double w)
{
    double gammln(double xx);
    return gammln(z)+gammln(w)-gammln(z+w);
}


// Pillai's helper
double C(int s, double m, double n)
{
//     // Components 1 and 2
//     double c1,c2,temp;
//     c1=c2=temp=0;

//     for (int i=1; i<s+1; i++)
//     {
//         c1 += gammln(0.5 * (2*m + 2*n + s + i + 2));
//         c2 += gammln(0.5 * (2*m + i + 1)) + gammln(0.5 * (2*n + i + 1)) + gammln(0.5 * i) ;
//     }
//     // ... and finally:
//     double out = log(pow(3.141593,0.5*s)) + c1 - c2;
//     return out;
}

// Function that will be integrated
double I(const double t, const double m, const double n)
{
    return pow(1/t,-1*m) * pow(1-t,n);
    //return pow(t,m) * pow(1-t,n);
    //return pow(t,2) * pow(1-t,991.5);
}

// Pillai's main function
long double pillai(int N, int p, int q, double lroot)
{

//     // Required input for Pillai's approximation. p must be < q; x is the
//     // largest root or eingenvalue (largest eigenvalue = largest cancor ^ 2)
//     int p2 = p <= q ? p : q;
//     int q2 = p > q ? p : q;
//     int s = p2;
//     double m = 0.5 * (q2-p2-1);
//     double n = 0.5 * (N-p2-q2-2) ;

//     double Csmn = C(s,m,n);
//     double Cs_1mn = C(s-1,m,n);

//     // Calculating log( V )
//     // We need 1 V for each of the 's-1' k's needed in the formula for P

//     vector_tld V(s-1);

//     for (int i=1; i<s; i++)
//     {
	
//         // i==1? If so, then the numerator of V = 1
//         if (i == 1)
//         {
//             double num = log(1);
//             double denom = Cs_1mn;
//             V[i-1] = num - denom;
// 	}
//         // for i>1
//         else
//         {
//             double a=1;
//             for (int j=1; j<i; j++)
//             {
//                 a *= (2*m+s-j+1) / (2*m+2*n+2*s-j+1);
//             }

// 	    double num = log( bico(s-1,i-1) * a );
// 	    double denom = Cs_1mn;

// 	    V[i-1] = num - denom ;
// 	}
//     }

//     // Computing k
//     vector_tld k(s,0);             // with 0s on it (the first element will NOT be used later on  
//     for (int i=1; i<s; i++)
//     {
//         k[i] = ( exp( Csmn + V[i-1] ) - (m+s-i+1)*k[i-1] )  /  (m+n+s-i+1);
//      }


//     // Computing the cdf
//     double sum_v;
//     //lroot=pow(0.05,2);  //test
//     long double cdf; 
//     for (int i=1; i<s; i++)
//     {
//  	sum_v += pow(-1.0,i) * k[i] * ( pow(lroot,m+s-i) * pow(1-lroot,n+1) );
//     }
    
//     // If min(p,q) is even...
//     if ( s%2 == 0 )
//     {
// 	cdf = 1 + sum_v;
//     }
//     // If s is odd
//      else
//     {
// 	double num = log( qromo(I,0.0,lroot,m,n) ); 
// 	double denom = betaln(m+1,n+1);
// 	double last = exp(num - denom);
// 	if (last > 1)
// 	    last = 1;
// 	cdf = last + sum_v;
//     }

//     return 1-cdf;

}


// Bartlet's test for all canonical correlations
long double bartlett(int N, int p, int q, vector_t eigen)
{
    int p2 = p <= q ? p : q; // Number of canonical correlations
    double prod_eigen=1;
    for (int i=0; i<p2; i++)
    {
	prod_eigen *= (1-eigen[i]);
    }
    double chisq = -1*(N - 1 - 0.5*(p+q+1)) * log(prod_eigen);
    double pvalue = chiprobP(chisq,p*q);
    return pvalue;
}

 


// Transpose of a matrix
void transposeMatrix(matrix_t & M)
{
    int rM = M.size();
    int cM = M[1].size();
    matrix_t tM;
    sizeMatrix(tM,cM,rM);
    for (int r=0; r<cM; r++)
    {
	for (int c=0; c<rM; c++)
	{
	    tM[r][c] = M[c][r];
	}
    }
    M = tM;
}


int calcGENEPIMeanVariance(vector<CSNP*> &,
			    int, int,
			    bool,
			    Plink *,
			    vector<double> &,
			    vector<vector<double> > &,
			    vector<Individual*>&,
			    vector<int> &,
			    vector<int> &);

void CCA_logit(bool perm, 
               vector<vector<int> > & blperm,
	       Set & S,
	       Plink & P);


void CCA_caseonly(bool perm,
		  vector<vector<int> > & blperm,
		  Set & S,
		  Plink & P);


void Plink::driverSCREEPI()
{

  ///////////////////////////////
  // Gene-based epistasis
  

  //////////////////////////////////////////
  // Case-control samples only

  affCoding(*this);


  //////////////////////////////////////////
  // SNP-major mode analysis

  if (!par::SNP_major)
    Ind2SNP();
  
  //////////////////////////////////////////
  // Requires that sets have been speciefied
  if (par::set_test) readSet();
  else error("Need to specify genes with --set {filename} when using --genepi\n");

    
  //////////////////
  // SET statistics

  Set S(snpset);


  //////////////////////////////////////////////
  // Prune SET (0-sized sets, MAF==0 SNPs, etc) 

  S.pruneSets(*this);
  
  int ns = snpset.size();
  if (ns < 2)
    error("Need to specify at least two fully valid sets\n");


  int n = 0;
  int ncase = 0;
  
  /////////////////////////////////////////////////////////
  // Prune based on VIF

  string original_outfile = par::output_file_name;

  // Case-control? Prune cases and controls together...
  if (!par::epi_caseonly)
  {   
      printLOG("\nConsidering cases and controls: ");
      setFlags(false);
      vector<Individual*>::iterator person = sample.begin();
      while ( person != sample.end() )
      {
	  if ( ! (*person)->missing )
	  {
	      (*person)->flag = true;
	      n++;
	  }
	  person++;
      }
  
      par::output_file_name += ".all";
      S.pruneMC(*this,false,par::vif_threshold);
      //S.pruneMC(*this,false,1000);
  }

  // Case-only? Prune cases only...
  else
  {
      printLOG("\nConsidering cases: ");
      setFlags(false);
      vector<Individual*>::iterator person = sample.begin();
      while ( person != sample.end() )
      {
	  if ( (*person)->aff && ! (*person)->missing )
	  {
	      (*person)->flag = true;
	      ncase++;
	  }
	  person++;
          n++;
      }

      par::output_file_name += ".case";
      S.pruneMC(*this,false,par::vif_threshold);
      //S.pruneMC(*this,false,1000);
  }

  par::output_file_name = original_outfile;

  // Write finalized set
  ofstream SET1, SET2;
  string f = par::output_file_name + ".all.set.in";
  
  printLOG("Writing combined pruned-in set file to [ " + f + " ]\n");
  SET1.open(f.c_str(),ios::out);

  f = par::output_file_name + ".all.set.out";
  printLOG("Writing combined pruned-out set file to [ " + f + " ]\n");
  SET2.open(f.c_str(),ios::out);

  for (int s=0; s<snpset.size(); s++)
  {
      
      int nss = snpset[s].size();
      
      SET1 << setname[s] << "\n";
      SET2 << setname[s] << "\n";
      
      for (int j=0; j<nss; j++)
      {
	  if (S.cur[s][j])
	      SET1 << locus[snpset[s][j]]->name << "\n";
	  else
	      SET2 << locus[snpset[s][j]]->name << "\n";
      }
      
      SET1 << "END\n\n";
      SET2 << "END\n\n";
  }
  
  SET1.close();
  SET2.close();
  

  // Prune empty sets once more:

  S.pruneSets(*this);
  
  ns = snpset.size();
  if (ns < 2)
      error("Need to specify at least two fully valid sets\n");


  ////////////////////////////////
  // Set up permutation structure

  // Specialized (i.e. cannot use Perm class) as this 
  // requires a block-locus permutation

  // First block is fixed
  
  vector<vector<int> > blperm(ns);
  vector<vector<int> > blperm_case(ns);
  vector<vector<int> > blperm_control(ns);

  for (int i=0; i<ns; i++)
  {
      // A slot for each individual per locus
      for (int j=0; j<n; j++)
	  if ( ! sample[j]->missing )
	      blperm[i].push_back(j);
 
      // A slot for each individual per locus
      for (int j=0; j<n; j++)
	  if ( ! sample[j]->missing && sample[j]->aff )
	      blperm_case[i].push_back(j);

      // A slot for each individual per locus
      for (int j=0; j<n; j++)
	  if ( ! sample[j]->missing && !sample[j]->aff )
	      blperm_control[i].push_back(j);
  }


  ////////////////////////////////////////////
  // Open file and print header for results

  ofstream EPI(f.c_str(), ios::out);
  EPI.open(f.c_str(), ios::out);
  EPI.precision(4);


  ////////////////////////////////////////
  // Analysis (calls genepi functions)

  if (!par::epi_caseonly)
      CCA_logit(false,blperm,S,*this);
  else
      CCA_caseonly(false,blperm_case,S,*this);

  if (!par::permute) 
   return;

  if (!par::silent)
    cout << "\n";


} // End of screepi



///////////////////////////
// CCA functions


///////////////////////////////////////////////////////////
// First CCA function: use for case-control logit analysis

void CCA_logit(bool perm, 
	       vector<vector<int> > & blperm,
	       Set & S,
	       Plink & P)
  
{

  ///////////////
  // Output results

      ofstream EPI;

      if (!perm)
      {  
	  string f = par::output_file_name+".genepi";
	  P.printLOG("\nWriting gene-based epistasis tests to [ " + f + " ]\n");
	  EPI.open(f.c_str(), ios::out);
	  EPI.precision(4);

	  EPI << setw(12) << "NIND"  << " "
	      << setw(12) << "GENE1"  << " "
	      << setw(12) << "GENE2"  << " "	 
	      << setw(12) << "NSNP1"  << " "
	      << setw(12) << "NSNP2"  << " "
	      << setw(12) << "P" << " "
	      << "\n";

      }


      //////////////////////////////////
      // Canonical correlation analysis

      int ns = P.snpset.size();

      // Consider each pair of genes
      for (int s1=0; s1 < ns-1; s1++)
      {
	for (int s2 = s1+1; s2 < ns; s2++)
        {


	    ////////////////////////////////////////////////////////
	    // Step 1. Construct covariance matrix (cases and controls together)
	    //    And partition covariance matrix:
	    //    S_11  S_21
	    //    S_12  S_22
	      
	    int n1=0, n2=0;
	      
	    vector<vector<double> > sigma(0);
	    vector<double> mean(0);
	    vector<CSNP*> pSNP(0);
	      
	    /////////////////////////////
	    // List of SNPs for both loci
	      
	    for (int l=0; l<P.snpset[s1].size(); l++)
            {
	      if ( S.cur[s1][l] )
	      {
		pSNP.push_back( P.SNP[ P.snpset[s1][l] ] );
		n1++;
	      }
	    }
	    for (int l=0; l<P.snpset[s2].size(); l++)
	    {		
              if ( S.cur[s2][l] )
	      {
		pSNP.push_back( P.SNP[ P.snpset[s2][l] ] );
		n2++;
              }
	    }

	    int n12 = n1 + n2;
	    int ne = n1 < n2 ? n1 : n2;
  
	    ///////////////////////////////////
	    // Construct covariance matrix (cases and controls together)
	      
	    P.setFlags(false);
	    vector<Individual*>::iterator person = P.sample.begin();

	    while ( person != P.sample.end() )
	      {
		(*person)->flag = true;
                person++;
	    }

	      
	    int nind = calcGENEPIMeanVariance(pSNP, 
					      n1,n2,
					      false,
					      &P,
					      mean, 
					      sigma, 
					      P.sample , 
					      blperm[s1],
					      blperm[s2] );
	    


	    ///////////////////////////
	    // Partition covariance matrix
	      
	    vector<vector<double> > I11;
	    vector<vector<double> > I11b;
	    vector<vector<double> > I12;
	    vector<vector<double> > I21;
	    vector<vector<double> > I22;
	    vector<vector<double> > I22b;
	      
	    sizeMatrix( I11, n1, n1);
            sizeMatrix( I11b, n1, n1);
	    sizeMatrix( I12, n1, n2);
	    sizeMatrix( I21, n2, n1);
	    sizeMatrix( I22, n2, n2);
	    sizeMatrix( I22b, n2, n2);             // For step 4b (eigenvectors for gene2)
	      
	    for (int i=0; i<n1; i++)
		for (int j=0; j<n1; j++)
                {
		  I11[i][j] = sigma[i][j];
		  I11b[i][j] = sigma[i][j];
                }
	      
	    for (int i=0; i<n1; i++)
		for (int j=0; j<n2; j++)
		    I12[i][j] = sigma[i][n1+j];
	      
	    for (int i=0; i<n2; i++)
		for (int j=0; j<n1; j++)
		    I21[i][j] = sigma[n1+i][j];
	      
	    for (int i=0; i<n2; i++)
		for (int j=0; j<n2; j++)
                {     
		  I22[i][j] = sigma[n1+i][n1+j];
	          I22b[i][j] = sigma[n1+i][n1+j];
		}


	    ////////////////////////////////////////////////////////
	    // Step 2. Calculate the p x p matrix M1 = inv(sqrt(sig11)) %*% sig12 %*% inv(sig22) %*% sig21 %*% inv(sqrt(sig11))
	    bool flag = true;
	    I11 = msqrt(I11);
	    I11 = svd_inverse(I11,flag);
       	    I22 = svd_inverse(I22,flag);

	    I22b = msqrt(I22b);// For Step 4b
	    I22b = svd_inverse(I22b,flag);
	    I11b = svd_inverse(I11b,flag);
      
	    matrix_t tmp;
	    matrix_t M1;
	      
	    multMatrix(I11, I12, tmp);
            multMatrix(tmp, I22, M1);
	    multMatrix(M1, I21, tmp);
	    multMatrix(tmp, I11, M1);


	    ////////////////////////////////////////////////////////
	    // Step 4a. Calculate the p eigenvalues and p x p eigenvectors of
	    // M (e). These are required to compute the coefficients used to
	    // build the p canonical variates a[k] for gene1 (see below)

            double max_cancor = 0.90;

            // Compute evalues and evectors
            Eigen gene1_eigen = eigenvectors(M1);

    
            // Sort evalues for gene 1. (the first p of these equal the first p of gene 2 (ie M2), if they are sorted)
            vector<double> sorted_eigenvalues_gene1 = gene1_eigen.d;
	    sort(sorted_eigenvalues_gene1.begin(),sorted_eigenvalues_gene1.end(),greater<double>());

            // Position of the largest canonical correlation that is <
            // max_cancor in the sorted vector of eigenvalues.  This will be
            // needed to use the right gene1 and gene2 coefficients to build
            // the appropriate canonical variates. 
            double cancor1=0;
            int cancor1_pos;          
            
            for (int i=0; i<n1; i++)
            {
              if ( sqrt(sorted_eigenvalues_gene1[i]) > cancor1 && sqrt(sorted_eigenvalues_gene1[i]) < max_cancor  )
              {
                cancor1 = sqrt(sorted_eigenvalues_gene1[i]);
                cancor1_pos = i;
                break;
              }
            }

            // Display largest canonical correlation and its position
	    //  cout << "Largest canonical correlation [position]\n"
	    //    << cancor1 << " [" << cancor1_pos << "]" << "\n\n" ;

            // Sort evectors. Rows must be ordered according to cancor value (highest first)
	    matrix_t sorted_eigenvectors_gene1 = gene1_eigen.z;
            vector<int> order_eigenvalues_gene1(n1);

	    for (int i=0; i<n1; i++)
            {
	      // Determine position of the vector associated with the ith cancor
              for (int j=0; j<n1; j++)
	      {
	        if (gene1_eigen.d[j]==sorted_eigenvalues_gene1[i])		 
                {
	          if (i==0)	   
                  {
		    order_eigenvalues_gene1[i]=j;
                    break;
                  }
                  else
                  {
		    if (j!=order_eigenvalues_gene1[i-1])
  		    {
                      order_eigenvalues_gene1[i]=j;
                      break;
                    }
	          }
                }
              }
	    }

            for (int i=0; i<n1; i++)
            {
		sorted_eigenvectors_gene1[i] = gene1_eigen.z[order_eigenvalues_gene1[i]];
	    }
	    //   cout << "Eigenvector matrix - unsorted:\n";
	    // display(gene1_eigen.z);
            //cout << "Eigenvector matrix - sorted:\n";
            //display(sorted_eigenvectors_gene1);


	    ////////////////////////////////////////////////////////
	    // Step 4b. Calculate the q x q eigenvectors of M2 (f). These are
	    // required to compute the coefficients used to build the p
	    // canonical variates b[k] for gene2 (see below). The first p are
	    // given by: f[k] = (1/sqrt(eigen[k])) * inv_sqrt_I22 %*% I21 %*%
	    // inv_sqrt_sig11 %*% e[k] for (k in 1:p) { e.vectors.gene2[,k] =
	    // (1/sqrt(e.values[k])) * inv.sqrt.sig22 %*% sig21 %*%
	    // inv.sqrt.sig11 %*% e.vectors.gene1[,k] }
           
             matrix_t M2;

             multMatrix(I22b, I21, tmp);
             multMatrix(tmp, I11b, M2);
             multMatrix(M2, I12, tmp);
             multMatrix(tmp, I22b, M2);
             Eigen gene2_eigen = eigenvectors(M2);
 
            //cout << "Eigenvalues Gene 2 - unsorted:\n";
            //display(gene2_eigen.d);
 
	    // Sort evalues for gene2
            vector<double> sorted_eigenvalues_gene2 = gene2_eigen.d;
            sort(sorted_eigenvalues_gene2.begin(),sorted_eigenvalues_gene2.end(),greater<double>());

            // Sort eigenvectors for gene2
            matrix_t sorted_eigenvectors_gene2 = gene2_eigen.z;
            vector<int> order_eigenvalues_gene2(n2);

            for (int i=0; i<n2; i++)
            {
		// Determine position of the vector associated with the ith cancor
		for (int j=0; j<n2; j++)
		{
		    if (gene2_eigen.d[j]==sorted_eigenvalues_gene2[i])
		    {
			if (i==0)
			{
			    order_eigenvalues_gene2[i]=j;
			    break;
			}
			else
			{
			    if (j!=order_eigenvalues_gene2[i-1])
			    {
				order_eigenvalues_gene2[i]=j;
				break;
			    }
			}
		    }
		}
            }

	    for (int i=0; i<n2; i++)
            {
                sorted_eigenvectors_gene2[i] = gene2_eigen.z[order_eigenvalues_gene2[i]];
            }

            //cout << "Eigenvector matrix Gene 2 - unsorted:\n";
	    //display(gene2_eigen.z);

            //cout << "Eigenvector matrix Gene 2 - sorted:\n";
            //display(sorted_eigenvectors_gene2);

            //exit(0);

            //////////////////////////////////////////////////////////////////////////////////
	    // Step 5 - Calculate the gene1 (pxp) and gene2 (pxq) coefficients
            // used to create the canonical variates associated with the p
            // canonical correlations

            transposeMatrix(gene1_eigen.z);
	    transposeMatrix(gene2_eigen.z);

	    matrix_t coeff_gene1;
            matrix_t coeff_gene2;

            multMatrix(gene1_eigen.z, I11, coeff_gene1);
	    multMatrix(gene2_eigen.z, I22b, coeff_gene2);

            //cout << "Coefficients for Gene 1:\n";
            //display(coeff_gene1);

            //cout << "Coefficients for Gene 2:\n";
            //display(coeff_gene2);

            //exit(0);

            ///////////////////////////////////////////////////////////////////////
            // Step 6 - Compute the gene1 and gene2 canonical variates
            // associated with the highest canonical correlation NOTE: the
            // original variables of data need to have the mean subtracted  first!
            // Otherwise, the resulting correlation between variate.gene1 and
            // variate.gene1 != estimated cancor.

            // For each individual, eg compos.gene1 =
            // evector.gene1[1]*SNP1.gene1 + evector.gene1[2]*SNP2.gene1 + ...

	    /////////////////////////////////
	    // Consider each SNP in gene1

	    vector<double> gene1(nind);
            
	    for (int j=0; j<n1; j++)
	    {

		CSNP * ps = pSNP[j];


           	///////////////////////////
		// Iterate over individuals

		for (int i=0; i< P.n ; i++)
		{

		    // Only need to look at one perm set
		    bool a1 = ps->one[i];
		    bool a2 = ps->two[i];

		    if ( a1 )
		    {
			if ( a2 ) // 11 homozygote
			{
			    gene1[i] += (1 - mean[j]) * coeff_gene1[order_eigenvalues_gene1[cancor1_pos]][j];
			}
                        else      // 12 
                        {
                            gene1[i] += (0 - mean[j]) * coeff_gene1[order_eigenvalues_gene1[cancor1_pos]][j];
                        }
		    }
		    else
		    {
                      if ( a2 )      // 21
                      {
                        gene1[i] += (0 - mean[j]) * coeff_gene1[order_eigenvalues_gene1[cancor1_pos]][j];
                      }
		      else           // 22 homozygote
		      {
			  gene1[i] += (-1 - mean[j]) * coeff_gene1[order_eigenvalues_gene1[cancor1_pos]][j];
		      }
		    }

		} // Next individual

	    } // Next SNP in gene1

           

            /////////////////////////////////
            // Consider each SNP in gene2
            vector<double> gene2(P.n);
            int cur_snp = -1;            
            for (int j=n1; j<n1+n2; j++)
            {

                cur_snp++;
                CSNP * ps = pSNP[j];


                
		// Iterate over individuals

		for (int i=0; i<P.n; i++)
		{
		    
		    // Only need to look at one perm set
		    bool a1 = ps->one[i];
		    bool a2 = ps->two[i];

		    if ( a1 )
		    {
			if ( a2 ) // 11 homozygote
			{
			    gene2[i] += (1 - mean[j]) * coeff_gene2[order_eigenvalues_gene2[cancor1_pos]][cur_snp];
			}
			else      // 12
			{
			    gene2[i] += (0 - mean[j]) * coeff_gene2[order_eigenvalues_gene2[cancor1_pos]][cur_snp];
			}
		    }
		    else
		    {
			if ( a2 )      // 21
			{
			    gene2[i] += (0 - mean[j]) * coeff_gene2[order_eigenvalues_gene2[cancor1_pos]][cur_snp];
			}
			else           // 22 homozygote
			{
			    gene2[i] += (-1 - mean[j]) * coeff_gene2[order_eigenvalues_gene2[cancor1_pos]][cur_snp];
			}
		    }
		    
		} // Next individual
		
	    } // Next SNP in gene2


            // Store gene1.variate and gene2.variate in the multiple_covariates field of P.sample
	    // TO DO: NEED TO CHECK IF FIELDS ARE EMPTY FIRST!

	    for (int i=0; i<P.n; i++)
	    {
		P.sample[i]->clist.resize(2);
		P.sample[i]->clist[0] = gene1[i];
		P.sample[i]->clist[1] = gene2[i];
	    }

            ///////////////////////////////////////////////
	    // STEP 7 - Logistic or linear regression epistasis test
	    //
	    
	    Model * lm;
	    	   
	    if (par::bt)
	    {
		LogisticModel * m = new LogisticModel(& P);
		lm = m;
	    }
	    else
	    {
		LinearModel * m = new LinearModel(& P);
		lm = m;
	    }

	    // No SNPs used
	    lm->hasSNPs(false);

	    // Set missing data
	    lm->setMissing();

 	    // Main effect of GENE1 1. Assumes that the variable is in position 0 of the clist vector
	    lm->addCovariate(0);
	    lm->label.push_back("GENE1");

	    // Main effect of GENE 2. Assumes that the variable is in position 1 of the clist vector
	    lm->addCovariate(1);
	    lm->label.push_back("GENE2");

	    // Epistasis
	    lm->addInteraction(1,2);
	    lm->label.push_back("EPI");

	    // Build design matrix
	    lm->buildDesignMatrix();

	    // Prune out any remaining missing individuals
// No longer needed (check)
//	    lm->pruneY();

	    // Fit linear model
	    lm->fitLM();


	    // Did model fit okay?
	    lm->validParameters();

	    // Obtain estimates and statistic
	    lm->testParameter = 3; // interaction
	    vector_t b = lm->getCoefs();
	    double chisq = lm->getStatistic();
	    double logit_pvalue = chiprobP(chisq,1);


	    // Clean up
	    delete lm;



            /////////////////////////////
            // OUTPUT

	    EPI << setw(12) << nind  << " "
		<< setw(12) << P.setname[s1]  << " "
                << setw(12) << P.setname[s2] << " "
                << setw(12) << n1  << " "
                << setw(12) << n2 << " "
                << setw(12) << logit_pvalue << " "
                << "\n";


        }  // End of loop over genes2
 

      }  // End of loop over genes1
      

      EPI.close();


}  // End of CCA_logit() 




///////////////////////////////////////////////////////////
// Second CCA function: use for case-control only

void CCA_caseonly(bool perm, 
	          vector<vector<int> > & blperm_case,
		  Set & S,
		  Plink & P)
  
{


    ///////////////
    // Output file

    ofstream EPI;

    if (!perm)
    {  
	string f = par::output_file_name+".genepi";
	P.printLOG("\nWriting gene-based epistasis tests to [ " + f + " ]\n");
	EPI.open(f.c_str(), ios::out);
	EPI.precision(4);

	EPI << setw(12) << "NIND"  << " "
	    << setw(12) << "GENE1"  << " "
	    << setw(12) << "GENE2"  << " "
	    << setw(12) << "NSNP1"  << " "
	    << setw(12) << "NSNP2"  << " "
	    << setw(12) << "CC1" << " "
	  //	    << setw(12) << "PILLAI" << " "
	    << setw(12) << "BART" << " "
	    << "\n";
    }

    //////////////////////////////////
    // Canonical correlation analysis

    // Number of genes
    int ns = P.snpset.size(); 


    // Consider each pair of genes
    for (int s1=0; s1 < ns-1; s1++)
    {
	for (int s2 = s1+1; s2 < ns; s2++)
        {


	    ////////////////////////////////////////////////////////
	    // Step 1. Construct covariance matrix (cases only)
	    //    And partition covariance matrix:
	    //    S_11  S_21
	    //    S_12  S_22
	          
	    int n1=0, n2=0;
	          
	    vector<vector<double> > sigma(0);
	    vector<double> mean(0);
	    vector<CSNP*> pSNP(0);


	    /////////////////////////////
	    // List of SNPs for both loci
	          
	    for (int l=0; l<P.snpset[s1].size(); l++)
            {
		if ( S.cur[s1][l] )
		{
		    pSNP.push_back( P.SNP[ P.snpset[s1][l] ] );
		    n1++;
		}
	    }
	    for (int l=0; l<P.snpset[s2].size(); l++)
	    {
		if ( S.cur[s2][l] )
		{
		    pSNP.push_back( P.SNP[ P.snpset[s2][l] ] );
		    n2++;
		}
	    }

	    // NOTE: we need to make sure that n1 < n2. Migth cause problems below if this is not the case.// *********
	    int n12 = n1 + n2;
	    int ne = n1 < n2 ? n1 : n2;// ne = min(p,q)

  
            ///////////////////////////////////////////////////
	    // Choose cases-only 
	          
	    P.setFlags(false);
	    vector<Individual*>::iterator person = P.sample.begin();
            int ncase=0;
	    while ( person != P.sample.end() )
	     {
		if ( (*person)->aff && !(*person)->missing ) 
		{
		    (*person)->flag = true;
		    ncase++;
		}
                person++;

	     }
	    
	    int nind = calcGENEPIMeanVariance(pSNP, 
					      n1,n2,
					      false,
					      &P,
					      mean, 
					      sigma, 
					      P.sample , 
					      blperm_case[s1],
					      blperm_case[s2] );

	    ///////////////////////////
	    // Partition covariance matrix
	          
	    vector<vector<double> > I11;
	    vector<vector<double> > I11b;
	    vector<vector<double> > I12;
	    vector<vector<double> > I21;
	    vector<vector<double> > I22;
	    vector<vector<double> > I22b;
	          
	    sizeMatrix( I11, n1, n1);
            sizeMatrix( I11b, n1, n1);
	    sizeMatrix( I12, n1, n2);
	    sizeMatrix( I21, n2, n1);
	    sizeMatrix( I22, n2, n2);
	    sizeMatrix( I22b, n2, n2);             // For step 4b (eigenvectors for gene2)
	          
	    for (int i=0; i<n1; i++)
		for (int j=0; j<n1; j++)
                {
		    I11[i][j] = sigma[i][j];
		    I11b[i][j] = sigma[i][j];
                }
	          
	    for (int i=0; i<n1; i++)
		for (int j=0; j<n2; j++)
		    I12[i][j] = sigma[i][n1+j];
	          
	    for (int i=0; i<n2; i++)
		for (int j=0; j<n1; j++)
		    I21[i][j] = sigma[n1+i][j];
	          
	    for (int i=0; i<n2; i++)
		for (int j=0; j<n2; j++)
                {     
		    I22[i][j] = sigma[n1+i][n1+j];
		    I22b[i][j] = sigma[n1+i][n1+j];
		}


	    ////////////////////////////////////////////////////////
	    // Step 2. Calculate the p x p matrix M1 = inv(sqrt(sig11)) %*% sig12 %*% inv(sig22) %*% sig21 %*% inv(sqrt(sig11))

	    bool flag = true;

	    I11 = msqrt(I11);
	    I11 = svd_inverse(I11,flag);
	    I22 = svd_inverse(I22,flag);

	    I22b = msqrt(I22b);// For Step 4b
	    I22b = svd_inverse(I22b,flag);
	    I11b = svd_inverse(I11b,flag);
      
	    matrix_t tmp;
	    matrix_t M1;
	          
	    multMatrix(I11, I12, tmp);
            multMatrix(tmp, I22, M1);
	    multMatrix(M1, I21, tmp);
	    multMatrix(tmp, I11, M1);


	    ////////////////////////////////////////////////////////
	    // Step 3. Determine the p eigenvalues of M1. The sqrt(eigen(M)) =
	    // p canonical correlations Identify the largest can corr 

            // Compute evalues and evectors
            vector_t eigen = eigenvalues(M1);

            // Sort eigenvalues
	    vector<double> sorted_eigen = eigen;
            sort(sorted_eigen.begin(),sorted_eigen.end(),greater<double>());

            // P-value
	    //	    long double pillai_pvalue = pillai(ncase,n1,n2,sorted_eigen[0]);
            long double bartlett_pvalue = bartlett(ncase,n1,n2,sorted_eigen);


            /////////////////////////////////////////////////////////////////////
	    // OUTPUT
	    
	    EPI << setw(12) << ncase  << " "
		<< setw(12) << P.setname[s1]  << " "
		<< setw(12) << P.setname[s2] << " "
		<< setw(12) << n1  << " "
		<< setw(12) << n2 << " "
                << setw(12) << sqrt(sorted_eigen[0]) << " "
	      //                << setw(12) << pillai_pvalue << " "
                << setw(12) << bartlett_pvalue << " "
		<< "\n";


	}  // End of loop over genes2

    }  // End of loop over genes1

    EPI.close();


}  // End of CCA_caseonly





//////////////////////////////////
// Helper functions

int calcGENEPIMeanVariance(vector<CSNP*> & pSNP,
			   int n1,
			   int n2,
			   bool perm, 
			   Plink * P,
			   vector<double> & mean,
			   vector<vector<double> > & variance,
			   vector<Individual*> & sample,
			   vector<int> & gp1,
			   vector<int> & gp2 )
  
{
  
  // Return number of individuals that the mean and variance matrix
  // are based on 

  bool casewise_deletion = false;


  // Calculate mean and variance for n1+n2 x n1+n2 matrix

  // Individual order in n1 , n2 deteremined by g1, g2
  // (i.e. block-based permutation)
  
  // Under permutations, mean and variances won't change
  // Store means only for now

  int nss = pSNP.size();
  
  // Original calculation?
  if (!perm)
    mean.resize(nss,0);
  
  vector<int> cnt(nss,0);
  variance.resize(nss);
  
  for (int j=0; j<nss; j++)
    variance[j].resize(nss,0);
  
      
  /////////
  // Means

  /////////////////////////////////
  // Consider each SNP in this set
  
  for (int j=0; j<nss; j++)
    {
      
      CSNP * ps = pSNP[j];
      
      ///////////////////////////
      // Iterate over individuals
      
      for (int i=0; i< P->n; i++)
	{	      

 	  // Only need to look at one perm set
	  bool a1 = ps->one[gp1[i]];
	  bool a2 = ps->two[gp2[i]];

	  if ( a1 )
	    {		   
	      if ( a2 ) // 11 homozygote
		{
		  mean[j]++;
		  cnt[j]++;
		}
	    }
	  else 
	    {
	      cnt[j]++;
	      if ( ! a2  ) // 00 homozygote
		mean[j]--;
	      
	    }
	  
	} // Next individual	  
      
    } // Next SNP in set
  
  for (int j=0; j<nss; j++)
    mean[j] /= (double)cnt[j];
  
  
  
  /////////////////////////////////////
  // Iterate over pairs of SNPs in SETs
  
  // First SNP 
  for (int j1=0; j1<nss; j1++)
    {
      CSNP * ps1 = pSNP[j1];

      // Second SNP
      for (int j2=0; j2<nss; j2++)
	{
	  CSNP * ps2 = pSNP[j2];

	  // Iterate over individuals
	  
	  for (int i=0; i<P->n; i++)
	    {
	      
	      bool a1, a2;
	      if (j1<n1)
		{
		  a1 = ps1->one[gp1[i]];
		  a2 = ps1->two[gp1[i]];
		}
	      else
		{
		  a1 = ps1->one[gp2[i]];
		  a2 = ps1->two[gp2[i]];
		}
	      
	      bool b1, b2;
	      if (j1<n1)
		{
		  b1 = ps2->one[gp1[i]];
		  b2 = ps2->two[gp1[i]];
		}
	      else
		{
		  b1 = ps2->one[gp2[i]];
		  b2 = ps2->two[gp2[i]];
		}
	      
	      
	      // Mean substitution
	      double v1=mean[j1], v2=mean[j2];
	      
	      // First SNP
	      if ( a1 )
		{		      
		  if ( a2 )    // 11 homozygote
		    {
		      v1 = 1;
		    }
		}
	      else 
		{
		  if ( ! a2  ) // 00 homozygote
		    {
		      v1 = -1;
		    }
		  else
		    v1 = 0;       // 01 heterozygote
		}
	      
	      
	      // Second SNP
	      if ( b1 )
		{		      
		  if ( b2 )    // 11 homozygote
		    {
		      v2 = 1;
		    }
		}
	      else 
		{
		  if ( ! b2  ) // 00 homozygote
		    {
		      v2 = -1;
		    }
		  else
		    v2 = 0;       // 01 heterozygote
		}
	      
	      
	      // Contribution to covariance term
	      variance[j1][j2] += ( v1 - mean[j1] ) * ( v2 - mean[j2] );
	      
		} // Next individual
	    } // Second SNP
	} // First SNP
        

      // Make symmetric covariance matrix
      for (int i=0; i<nss; i++)
	for (int j=i; j<nss; j++)
	  {
	    variance[i][j] /= (double)(P->n);
 	    variance[j][i] = variance[i][j];
 	  }
      

      // Mean-imputation uses everybody

      return P->n;
      
 }
  
 

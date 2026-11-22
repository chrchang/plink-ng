

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
#include <cmath>

#include "plink.h"
#include "perm.h"
#include "options.h"
#include "helper.h"
#include "stats.h"


void Plink::calcMH()
{

  ///////////////////////////////////
  // Basic 2 x 2 x K CMH test
  // i.e. Disease x allele x strata
  // is taken care of in assoc.cpp 
  // (i.e. allows for permutation, sets, etc)


  if (!par::SNP_major) Ind2SNP();


  //////////////////////////////////
  // Any individual not assigned to a cluster, 
  // making missing phenotype
  
  vector<Individual*>::iterator person = sample.begin();
  while ( person != sample.end() )
    {
      if ( (*person)->sol < 0 ) 
	(*person)->missing = true;
      person++;
    }


  ///////////////////////////////////
  // Generalized I x J x K CMH test

  // Either ordinal or normal
  // i.e. test strata X SNP controlling for disease
  
  if (par::CMH_test_2 || par::CMH_test_ORD )
    {

      if (par::CMH_test_ORD && !par::bt) 
      error("--mh-ord specified but the phenotype is only binary: use --mh");

      if (nk==1) 
	error("No clusters defined for --mh2 test, i.e. K=1");

      string f = par::output_file_name + ".cmh2";
	if (par::CMH_test_ORD)
        f = par::output_file_name + ".cmh.ord";

      ofstream MHOUT;
      MHOUT.open(f.c_str(),ios::out);
      
      MHOUT << setw(4) << "CHR" << " "
	      << setw(par::pp_maxsnp) << "SNP" << " "
	      << setw(10) << "CHISQ" << " "
	      << setw(10) << "P" << "\n";
	        
      MHOUT.precision(4);
      
      if (par::CMH_test_ORD)
	{
	  printLOG("Cochran-Mantel-Haenszel IxJxK ordinal test, K = " 
		   + int2str(nk) + "\n");
	  printLOG("Testing SNP x ORDINAL DISEASE | STRATUM (option --mh-ord)\n");
	} 
      else 
	{ 
	  printLOG("Cochran-Mantel-Haenszel IxJxK test, K = " 
		   + int2str(nk) + "\n");
	  printLOG("Testing SNP x STRATUM | DISEASE (option --mh2)\n");
	}
      printLOG("Writing results to [ " + f + " ]\n");
      
      vector<CSNP*>::iterator s = SNP.begin();
      int l=0;
      while ( s != SNP.end() )
	{
	  
	  ///////////////////////// 	      
	  // Autosomal or haploid?
	  
	  bool Xchr=false, haploid=false;
	  if (par::chr_sex[locus[l]->chr]) Xchr=true;
	  else if (par::chr_haploid[locus[l]->chr]) haploid=true;
	  
	  if (haploid || Xchr ) 
	    error("--mh2 / --mh-ord cannot handle X/Y markers currently...");
	  
	  vector<int> X(0);  // SNP
	  vector<int> Y(0);  // Cluster
	  vector<int> Z(0);  // Phenotype
	  
	  vector<Individual*>::iterator person = sample.begin();
	  vector<bool>::iterator i1 = (*s)->one.begin();
	  vector<bool>::iterator i2 = (*s)->two.begin();
	  
	  while ( person != sample.end() )
	    {
	      
	      if ((*person)->missing) 
		{
		  // Next person
		  person++;
		  i1++;
		  i2++;
		  continue;
		}            	      
	      
	      // Only consider individuals who have been assigned to a cluster
	      if ( (*person)->sol >= 0 )
		{
		  if ( (!(*i1)) && (!(*i2)) ) 
		    { 
		      X.push_back(1); 
		      X.push_back(1); 
		    }
		 else if ( (!(*i1)) && *i2 ) 
		   { 
		     X.push_back(1); 
		     X.push_back(2); 
		   }
		  else if ( *i1 && *i2 ) 
		    { 
		      X.push_back(2); 
		      X.push_back(2); 
		    }
		  else 
		    {
		      // Next person
		      person++;
		      i1++;
		      i2++;
		      continue;
		    }
		  
		  Y.push_back((*person)->sol);
		  Y.push_back((*person)->sol);
		  
		  if (par::CMH_test_ORD) 
		    Z.push_back( (int)(*person)->phenotype );
		  else 
		    {
		      if ((*person)->phenotype==2) 
			{ 
			  Z.push_back(2); 
			  Z.push_back(2);
			} 
		      else { 
			Z.push_back(1); 
			Z.push_back(1); 
		      }       
		    }
		}
	      
	      // Next person
	      person++;
	      i1++;
	      i2++;
	      
	    }
	  
	  
	  vector<double> res;
	 
	  if ( par::CMH_test_ORD ) 
	    res = calcMantelHaenszel_ORD(X,Z,Y);
	  else 
	    res = calcMantelHaenszel_IxJxK(X,Y,Z);
	  
	  MHOUT << setw(4) << locus[l]->chr << " "
		<< setw(par::pp_maxsnp) << locus[l]->name << " "
		<< setw(10) << res[0] << " "
		<< setw(10) << chiprobP(res[0],res[1]) << "\n";
	  
	  // Next SNP
	  s++;
	  l++;
	}
      
      MHOUT.close();
    }
  
}


vector<double> Plink::calcMantelHaenszel_2x2xK(Perm & perm, bool original)
{

  // Should we perform BD test (K>1)
  if (nk<2) par::breslowday = false;

  ofstream MHOUT;

  if ( original )
    {
      
      //////////////////////////////////
      // Any individual not assigned to a cluster, making missing
      // phenotype (only need to do this once, for original)
      
      vector<Individual*>::iterator person = sample.begin();
      while ( person != sample.end() )
	{
	  if ( (*person)->sol < 0 ) 
	    (*person)->missing = true;
	  person++;
	}
      
      
      string f = par::output_file_name + ".cmh";
      MHOUT.open(f.c_str(),ios::out);
      
      MHOUT << setw(4) << "CHR" << " "
	    << setw(par::pp_maxsnp) << "SNP" << " "
	    << setw(10) << "BP" << " "
	    << setw(4) << "A1" << " " 
	    << setw(8) << "MAF" << " "
	    << setw(4) << "A2" << " "	
	    << setw(10) << "CHISQ" << " "
	    << setw(10) << "P" << " "
	    << setw(10) << "OR" << " "
	    << setw(10) << "SE" << " "
	    << setw(10) << string("L"+dbl2str(par::ci_level*100)) << " "  
	    << setw(10) << string("U"+dbl2str(par::ci_level*100)) << " ";    
      
      
      if (par::breslowday)
	MHOUT << setw(10) << "CHISQ_BD" << " "
	      << setw(10) << "P_BD" << " ";
      
      MHOUT << "\n";
      
      MHOUT.precision(4);
      
      printLOG("Cochran-Mantel-Haenszel 2x2xK test, K = " + int2str( nk) + "\n");

      if (par::breslowday)
	printLOG("Performing Breslow-Day test of homogeneous odds ratios\n");
      
      printLOG("Writing results to [ " + f + " ]\n");
      
      // Warnings, 
      if (par::breslowday && nk>10) 
	printLOG("** Warning ** Breslow-Day statistics require large N per cluster ** \n");
      
    }


  double zt = ltqnorm( 1 - (1 - par::ci_level) / 2  ) ; 
  

  // Cochran-Mantel-Haenszel 2x2xK test
  
  vector<double> results(nl_all);

  vector<CSNP*>::iterator s = SNP.begin();
  int l=0;
  while ( s != SNP.end() )
    {
      
      // Skip possibly
      if (par::adaptive_perm && !perm.snp_test[l])
	{
	  s++;
	  l++;
	  continue;
	}

      // Disease X allele X strata
      
      // Calculate mean of 11 cell for each strata
      vector<double> mean_11(nk,0);
      vector<double> var_11(nk,0);
 
      // Calculate statistic
      vector<double> n_11(nk,0);
      vector<double> n_12(nk,0);
      vector<double> n_21(nk,0);
      vector<double> n_22(nk,0);

      // Disease marginals
      vector<double> n_1X(nk,0); // disease
      vector<double> n_2X(nk,0); // no disease

      vector<double> n_X1(nk,0); // F allele
      vector<double> n_X2(nk,0); // T allele
      
      vector<double> n_TT(nk,0); // Total allele count
      

      ///////////////////////// 	      
      // Autosomal or haploid?
	      
      bool X=false, haploid=false;
      if (par::chr_sex[locus[l]->chr]) X=true;
      else if (par::chr_haploid[locus[l]->chr]) haploid=true;

      ////////////////////////
      // Consider each person
      
      vector<bool>::iterator i1 = (*s)->one.begin();
      vector<bool>::iterator i2 = (*s)->two.begin();
      vector<Individual*>::iterator gperson = sample.begin();

      while ( gperson != sample.end() )
	{
	  Individual * pperson = (*gperson)->pperson;

	  bool s1 = *i1;
	  bool s2 = *i2;
	  
	  // Affected individuals
	  if ( pperson->aff && !pperson->missing )
	    {

		// Haploid?
		if ( haploid || ( X && (*gperson)->sex ) ) 
		{

	      // Allelic marginal
	      if ( ! s1  )
		{
		  // FF hom
		      n_11[ pperson->sol ] ++ ;
		      n_X1[ pperson->sol ] ++ ;		    
		}
	      else 
		{
		  if ( ! s2 ) // FT
		    {
		      gperson++;
		      i1++;
		      i2++;
		      continue;  // skip missing genotypes
		    }
		  else // TT
		    {
		      n_12[ pperson->sol ] ++ ;
		      n_X2[ pperson->sol ] ++ ;		      
		    }
		}

	      // Disease marginal
	      n_1X[ pperson->sol ] ++;
	      n_TT[ pperson->sol ] ++;

		}
		else   // autosomal
		{
		
	      // Allelic marginal
	      if ( ! s1  )
		{
		  if ( ! s2 ) // FF hom
		    {
		      n_11[ pperson->sol ] +=2 ;
		      n_X1[ pperson->sol ] +=2 ;
		    }
		  else
		    {
		      n_11[ pperson->sol ]++ ; // FT het
		      n_12[ pperson->sol ]++ ;
		      n_X1[ pperson->sol ]++ ;
		      n_X2[ pperson->sol ]++ ;		      
		    }
		}
	      else 
		{
		  if ( ! s2 ) // FT
		    {
		      gperson++;
		      i1++;
		      i2++;
		      continue;  // skip missing genotypes
		    }
		  else // TT
		    {
		      n_12[ pperson->sol ] +=2 ;
		      n_X2[ pperson->sol ] +=2 ;		      
		    }
		}

	      // Disease marginal
	      n_1X[ pperson->sol ] += 2;
	      n_TT[ pperson->sol ] += 2;

	      } // end autosomal

	    }
	  else if ( ! pperson->missing ) // Unaffecteds
	    {

		// Haploid?
		if ( haploid || ( X && (*gperson)->sex ) ) 
		{

	      // Allelic marginal
	      if ( ! s1  )
		{
		  // FF hom
		      n_21[ pperson->sol ] ++ ;
		      n_X1[ pperson->sol ] ++ ;		    
		}
	      else 
		{
		  if ( ! s2 ) // FT
		    {
		      gperson++;
		      i1++;
		      i2++;
		      continue;  // skip missing genotypes
		    }
		  else // TT
		    {
		      n_22[ pperson->sol ] ++ ;
		      n_X2[ pperson->sol ] ++ ;		      
		    }
		}

	      // Disease marginal
	      n_2X[ pperson->sol ] ++;
	      n_TT[ pperson->sol ] ++;

		}
		else   // autosomal 
		{
	      // Allelic marginal
	      if ( ! s1 )
		{
		  if ( ! s2 ) // FF
		    {
		      n_X1[ pperson->sol ] +=2 ;
		      n_21[ pperson->sol ] +=2 ;
		    }
		  else
		    {
		      n_X1[ pperson->sol ] ++ ;
		      n_X2[ pperson->sol ] ++ ;		      
		      n_21[ pperson->sol ] ++ ;
		      n_22[ pperson->sol ] ++ ;
		    }
		}
	      else 
		{
		  if ( ! s2 ) // FT
		    {
		      gperson++;
		      i1++;
		      i2++;
		      continue;  // skip missing genotypes
		    }
		  else // TT
		    {
		      n_X2[ pperson->sol ] +=2 ;		      		      
		      n_22[ pperson->sol ] +=2 ;
		    }
		}     

	      // disease marginal
	      n_2X[ pperson->sol ] += 2;
	      n_TT[ pperson->sol ] += 2;

		}  // end autosomal			      
	    } // end unaffected

	  gperson++;
	  i1++;
	  i2++;
	  
	} // count next individual
    
      
      
      

      // Finished iterating over individuals: cluster needs at least 2 
      // nonmissing individuals
      
      vector<bool> validK(nk,false);
      for (int k=0; k<nk; k++)
	if (n_TT[k]>=2) validK[k]=true;
      
      for (int k=0; k<nk; k++)
	{
	  if (validK[k])
	    {
	      mean_11[k] = ( n_X1[k] * n_1X[k] ) / n_TT[k] ;
	      var_11[k] = ( n_X1[k] * n_X2[k] * n_1X[k] * n_2X[k] ) 
		/ ( n_TT[k]*n_TT[k]*(n_TT[k]-1) );

// 	      cout << k << " " 
// 		   << n_11[k] << " " 
// 		   << n_12[k] << " " 
// 		   << n_21[k] << " " 
// 		   << n_22[k] << "\n";
	      
	    }
	}

      
      double CMH = 0;
      double denom = 0;
      for (int k=0; k<nk; k++)
	{
	  if (validK[k])
	    {
	      CMH += n_11[k] - mean_11[k];
	      denom += var_11[k];
	    }
	}

      CMH *= CMH;
      CMH /= denom;
 
      // MH Odds ratio & CI
      double R = 0;
      double S = 0;
      vector<double> r2(nk);
      vector<double> s2(nk);      

      for (int k=0; k<nk; k++)
	{
	  if (validK[k])
	    {
	      r2[k] = (n_11[k]*n_22[k]) / n_TT[k];
	      s2[k] = (n_12[k]*n_21[k]) / n_TT[k];
	      R += r2[k]; 
	      S += s2[k];
	    }
	}
      double OR = R / S ;

      double v1 = 0, v2 = 0, v3 = 0;
      for (int k=0; k<nk; k++)
	{
	  if (validK[k])
	    {
	      v1 += (1/n_TT[k]) * ( n_11[k] + n_22[k] ) * r2[k] ;
	      v2 += (1/n_TT[k]) * ( n_12[k] + n_21[k] ) * s2[k] ;
	      v3 += (1/n_TT[k]) * ( ( n_11[k] + n_22[k] ) * s2[k] 
				+ ( n_12[k] + n_21[k] ) * r2[k] );
	    }
	}
      
      double SE = ( 1/(2*R*R) )  * v1
	+ (1/(2*S*S)) * v2
	+ (1/(2*R*S)) * v3 ;
      SE = sqrt(SE);
      
      double OR_lower = exp( log(OR) - zt * SE );
      double OR_upper = exp( log(OR) + zt * SE );
      
      if ( original )
	{

	  double pvalue = chiprobP(CMH,1);
	  
	  // Skip?, if filtering p-values
	  if ( par::pfilter && ( pvalue > par::pfvalue || pvalue < 0 ) ) 
	    goto skip_p_cmh;
	  

	  MHOUT << setw(4) << locus[l]->chr << " "
		<< setw(par::pp_maxsnp) << locus[l]->name << " "
		<< setw(10) << locus[l]->bp << " "
		<< setw(4) << locus[l]->allele1 << " " 
		<< setw(8) << locus[l]->freq << " " 
		<< setw(4) << locus[l]->allele2 << " ";
	  

	  if (realnum(CMH))
	    MHOUT << setw(10) << CMH << " "
		  << setw(10) << chiprobP(CMH,1) << " ";
	  else
	    MHOUT << setw(10) << "NA" << " "
		  << setw(10) << "NA" << " ";
	  
	  
	  if (realnum(OR))
	    MHOUT << setw(10) << OR << " ";
	  else
	    MHOUT << setw(10) << "NA" << " ";

	  if (realnum(SE))
	    MHOUT << setw(10) << SE << " ";
	  else 
	    MHOUT << setw(10) << "NA" << " ";

	  if (realnum(OR_lower))
	    MHOUT << setw(10) << OR_lower << " ";
	  else
	    MHOUT << setw(10) << "NA" << " ";

	  if (realnum(OR_upper))
	    MHOUT << setw(10) << OR_upper << " ";
	  else
	    MHOUT << setw(10) << "NA" << " ";
	  

	  // Optional Breslow-Day test of homogeneity of odds ratios
	  if (par::breslowday)
	    {
	      
	      double amax;
	      double bb;
	      double determ;
	      double as_plus;
	      double as_minus;
	      double Astar;
	      double Bstar;
	      double Cstar;
	      double Dstar;
	      double Var;
	      double BDX2 = 0;
	      int df = 0;
	      for (int k=0; k<nk; k++)
		{
		  if (validK[k])
		    {
		      df++;
		      amax = (n_1X[k] < n_X1[k]) ? n_1X[k] : n_X1[k];
		      bb = n_2X[k] + n_1X[k] * OR - n_X1[k] * (1-OR);
		      determ = sqrt(bb*bb + 4*(1-OR) * OR * n_1X[k] * n_X1[k]);
		      as_plus = ( -bb + determ ) / ( 2 - 2 * OR );
		      as_minus = ( -bb - determ ) / ( 2 - 2 * OR );		      
		      Astar =  as_minus <= amax && as_minus >= 0 ? as_minus  : as_plus ;
		      Bstar = n_1X[k] - Astar;
		      Cstar = n_X1[k] - Astar;
		      Dstar = n_2X[k] - n_X1[k] + Astar;
		      Var = 1/(1/Astar + 1/Bstar + 1/Cstar + 1/Dstar);
		      BDX2 += ( (n_11[k] - Astar) * ( n_11[k] - Astar ) ) / Var ; 
		    }
		}
	      
	      double BDp = chiprobP( BDX2 , df-1 ); 
	      
	      if ( BDp > -1 )   
		MHOUT << setw(10) << BDX2 << " " 
		      << setw(10) << BDp << " ";	  
	      else
		MHOUT << setw(10) << "NA" << " " 
		      << setw(10) << "NA" << " ";	  
	      
	    }
	  
	  MHOUT << "\n";
	  
	}

    skip_p_cmh:
      
      // Store for permutation procedure, based 2x2xK CMH result
      results[l] = CMH;
      
      // Next SNP
      s++;
      l++;

    }
  

  if (original)
    MHOUT.close();

  return results;
  
}

  
vector<double> Plink::calcMantelHaenszel_IxJxK(vector<int> & X,
					       vector<int> & Y,
					       vector<int> & Z)
{
  
  if (X.size() != Y.size() || Y.size() != Z.size() || X.size() != Z.size() )
    error("Internal problem:\n  problem in calcMantelHaenszel_IxJxK()...uneven input columns");
  
  // Determine unique elements
  int nx=0, ny=0, nz=0;
  map<int,int> mx;
  map<int,int> my;
  map<int,int> mz;

  for (unsigned int i=0; i<X.size(); i++)
    {
      if (mx.find(X[i]) == mx.end())
	mx.insert(make_pair(X[i],nx++));

      if (my.find(Y[i]) == my.end())
	my.insert(make_pair(Y[i],ny++));

      if (mz.find(Z[i]) == mz.end())
	mz.insert(make_pair(Z[i],nz++));
    }

  // Generic function to calculate generalized IxJxK CMH
  // Assumes no missing data

  vector< vector<double> > N(nz);             // observed counts
  vector< vector<double> > U(nz);             // expected
  vector< vector< vector<double> > > V(nz);   // variance matrix
  
  vector<vector<double> > Tx(nz); // marginal totals
  vector<vector<double> > Ty(nz); // ..
  vector<double> T(nz);           // totals (per K)
  
  for (int k=0; k<nz; k++)
    {
      Tx[k].resize(nx);
      Ty[k].resize(ny);
      
      N[k].resize((nx-1)*(ny-1));
      U[k].resize((nx-1)*(ny-1));
      V[k].resize((nx-1)*(ny-1));
      for (int k2=0; k2<(nx-1)*(ny-1); k2++)
	{
	  N[k][k2] = U[k][k2] = 0;
	  V[k][k2].resize((nx-1)*(ny-1));
	  for (int k3=0; k3<(nx-1)*(ny-1); k3++)
	    V[k][k2][k3] = 0;
	}
    }
  

  // Consider each observation
  for (int i=0; i<X.size(); i++)
    {
      int vx = mx.find(X[i])->second;
      int vy = my.find(Y[i])->second;
      int vz = mz.find(Z[i])->second;
      
      // exclude nx + ny (upper limits)
      if (vx<nx-1 && vy<ny-1)
	N[vz][ vx + vy*(nx-1) ]++;
      
      Tx[vz][vx]++;
      Ty[vz][vy]++;
      T[vz]++;
    }
  
  
  // Determine valid clusters (at least 2 people)
  vector<bool> validK(nk,false);
  for (int k=0; k<nk; k++)
    if (T[k]>=2) validK[k]=true;
  

  // Calculate expecteds
  for (int k=0; k<nz; k++)
    {
      if (validK[k])
	{
	  for (int ix=0; ix<nx-1; ix++)
	    for (int iy=0; iy<ny-1; iy++)
	      {
		U[k][ix+iy*(nx-1)] = ( Tx[k][ix] * Ty[k][iy] ) / T[k];
		
		for (int ix2=0; ix2<nx-1; ix2++)
		  for (int iy2=0; iy2<ny-1; iy2++)
		    {
		      int dx=0;
		      int dy=0;
		      if (ix==ix2) dx=1;
		      if (iy==iy2) dy=1;
		      V[k][ix + iy*(nx-1)][ix2 + iy2*(nx-1)] = ( ( Tx[k][ix] * ( dx * T[k] - Tx[k][ix2] ) 
								   * Ty[k][iy] * ( dy *T[k] - Ty[k][iy2] )  )
								 / ( T[k]*T[k]*(T[k]-1) ) );
		      if (ix==ix2 && iy==iy2) 
			V[k][ix + iy*(nx-1)][ix2 + iy2*(nx-1)] 
			  = abs(V[k][ix + iy*(nx-1)][ix2 + iy2*(nx-1)]);
		    }
	      }
	}
    }
  

  vector<vector<double> > V0((nx-1)*(ny-1));
  for (int k2=0; k2<(nx-1)*(ny-1); k2++)
    V0[k2].resize((nx-1)*(ny-1));
  vector<double> N0((nx-1)*(ny-1));
  vector<double> U0((nx-1)*(ny-1));

  // Sum N, U and V over K
  for (int k=0; k<nz; k++)
    {
      if (validK[k])
	{
	  for (int i=0; i<(nx-1)*(ny-1); i++)
	    {
	      N0[i] += N[k][i];
	      U0[i] += U[k][i];
	      
	      for (int i2=0; i2<(nx-1)*(ny-1); i2++)
		V0[i][i2] += V[k][i][i2];
	    }
	}
    }
  
  
  bool flag = true;
  vector<double> tmp1((nx-1)*(ny-1),0);
  vector<double> tmp2((nx-1)*(ny-1),0);
  V0 = svd_inverse(V0,flag);
  for (int i=0; i<(nx-1)*(ny-1); i++)
    tmp1[i] = N0[i] - U0[i];
  
  // Matrix mult -- rows by columns
  
  for (int i=0; i<(nx-1)*(ny-1); i++)
    for (int j=0; j<(nx-1)*(ny-1); j++)
      tmp2[j] += tmp1[i] * V0[i][j];
 

  vector<double> result(2);

  // CMH Chi-square
  result[0]=0;
  for (int i=0; i<(nx-1)*(ny-1); i++)
    result[0] += tmp2[i] * tmp1[i];
  
  // DF
  result[1] = (nx-1)*(ny-1);
  return result;

}




void Plink::calcHomog()
{
  
  if (!par::SNP_major) Ind2SNP();

  string f = par::output_file_name + ".homog";
  ofstream MHOUT;
  MHOUT.open(f.c_str(),ios::out);

  MHOUT.precision(4);

  if (nk==0) error("No clusters (K=0)... cannot perform CMH tests");
 
  printLOG("Homogeneity of odds ratio test, K = " + int2str(nk) + "\n");

  if (nk<2) 
    {
      printLOG("** Warning ** less then 2 clusters specified... \n");
      printLOG("              cannot compute between-cluster effects ** \n");
      return;
    }

  if (nk>10) 
    printLOG("** Warning ** statistics can be unreliable if strata have small N ** \n");
  

  printLOG("Writing results to [ " + f + " ]\n");
  
  

  MHOUT << setw(4) << "CHR" << " " 
	<< setw(par::pp_maxsnp) << "SNP" << " " 
	<< setw(4) << "A1" << " "
	<< setw(4) << "A2" << " "
	<< setw(8) << "F_A" << " "
	<< setw(8) << "F_U" << " "
	<< setw(8) << "N_A" << " "
	<< setw(8) << "N_U" << " "
    	<< setw(8) << "TEST" << " "
	<< setw(10) << "CHISQ" << " "
	<< setw(4)  << "DF" << " "
	<< setw(10) << "P" << " "
	<< setw(10) << "OR" << "\n";
  
  
  ///////////////////////////////////
  // Create boolean affection coding

  affCoding(*this);


  //////////////////////////////////
  // Any individual not assigned to a cluster, 
  // making missing phenotype
  
  vector<Individual*>::iterator person = sample.begin();
  while ( person != sample.end() )
    {
      if ( (*person)->sol < 0 ) 
	(*person)->missing = true;
      person++;
    }

  	
  ///////////////////////////////
  // Iterate over SNPs
  
  vector<CSNP*>::iterator s = SNP.begin();
  int l=0;
  
  while ( s != SNP.end() )
    {	
      
      // Uncomment this if we allow permutation for the CMH
      // tests

      // In adaptive mode, possibly skip this test
      //      if (par::adaptive_perm && (!perm.snp_test[l])) 
      //	{
      //	  s++;
      //	  l++;
      //	  continue;
      //	}
      
      // Calculate statistic
      vector<double> n_11(nk,0);
      vector<double> n_12(nk,0);
      vector<double> n_21(nk,0);
      vector<double> n_22(nk,0);

      vector<double> lnOR(nk,0);
      vector<double> SEsq(nk,0);

      ///////////////// 	      
      // Autosomal or haploid?
      
      bool X=false, haploid=false;
      if (par::chr_sex[locus[l]->chr]) X=true;
      else if (par::chr_haploid[locus[l]->chr]) haploid=true;

      /////////////////////////////
      // Iterate over individuals
      
      vector<bool>::iterator i1 = (*s)->one.begin();
      vector<bool>::iterator i2 = (*s)->two.begin();
      vector<Individual*>::iterator gperson = sample.begin();
      

      while ( gperson != sample.end() )
	{
	  
	  // Phenotype for this person (i.e. might be permuted)
	  Individual * pperson = (*gperson)->pperson;
	  
	  // SNP alleles
	  
	  bool s1 = *i1;
	  bool s2 = *i2;

	  int hom = 2;
	  if ( haploid || ( X && (*gperson)->sex ) )
	    hom = 1;
	  
	  // Affected individuals
	  if ( pperson->aff && !pperson->missing )
	    {
	      
	      // Allelic marginal
	      if ( !s1 )
		{
		  if ( !s2 ) // FF hom
		    {
		      n_11[ pperson->sol ] += hom ;
		      
		    }
		  else
		    {
		      n_11[ pperson->sol ]++ ; // FT het
		      n_12[ pperson->sol ]++ ;
		    }
		}
	      else 
		{
		  if ( !s2 ) // FT
		    {
		      gperson++;
		      i1++;
		      i2++;
		      continue;  // skip missing genotypes
		    }
		  else // TT
		    {
		      n_12[ pperson->sol ] += hom ;
		    }
		}
	    }
	  else if ( ! pperson->missing ) // Unaffecteds
	    {
	      // Allelic marginal
	      if ( ! s1 )
		{
		  if ( ! s2 ) // FF
		    {
		      n_21[ pperson->sol ] += hom ;
		    }
	      else
		{
		  n_21[ pperson->sol ] ++ ;
		  n_22[ pperson->sol ] ++ ;
		}
		}
	      else 
		{
		  if ( ! s2 ) // FT
		    {
		      gperson++;
		      i1++;
		      i2++;
		      continue;  // skip missing genotypes
		    }
		  else // TT
		    {
		      n_22[ pperson->sol ] += hom ;
		    }
		}     
	    }
	 
	  // Next individual
	  gperson++;
	  i1++;
	  i2++;
 
	} 
      
      
      // Calculate log(OR) and SE(ln(OR)) for eacsh strata

      double X_total = 0;
      double X_assoc1 = 0;
      double X_assoc2 = 0;
      vector<double> X_indiv(nk,0);
     
      for (int k=0; k<nk; k++)
	{
	  
	  // Add 0.5 to each cell to reduce bias

  	  n_11[k] += 0.5;
  	  n_12[k] += 0.5;
  	  n_21[k] += 0.5;
  	  n_22[k] += 0.5;

	  // ln(OR)
	  
	  lnOR[k] = log ( ( n_11[k] * n_22[k] ) / ( n_12[k] * n_21[k] ) );
	  SEsq[k] =  1/n_11[k] + 1/n_12[k] + 1/n_21[k] + 1/n_22[k] ;
	  
	  X_indiv[k] = (lnOR[k] * lnOR[k]) / SEsq[k];
	  X_total += X_indiv[k];

	  // For the common, strata-adjusted test
	  X_assoc1 += lnOR[k] / SEsq[k];
	  X_assoc2 += 1/ SEsq[k];
	}
      
      // X_total is total chi-square on nk df
      // X_indiv are individual chi-squares, each on 1 df
      // X_homog is test for homogeneity of OR, with nk-1 df 
      // X_assoc is strata-adjusted test, with 1 df

      double X_assoc = (X_assoc1*X_assoc1)/X_assoc2;
      double X_homog = X_total - X_assoc;

      MHOUT << setw(4) << locus[l]->chr << " " 
	    << setw(par::pp_maxsnp) << locus[l]->name << " "
	    << setw(4) << locus[l]->allele1 << " "
	    << setw(4) << locus[l]->allele2 << " "
	    << setw(8) << "NA" << " "
	    << setw(8) << "NA" << " "
	    << setw(8) << "NA" << " "
	    << setw(8) << "NA" << " "	
	    << setw(6) << "TOTAL" << " "
	    << setw(10) << X_total << " "
	    << setw(4) << nk << " "
	    << setw(10) << chiprobP(X_total,nk) << " "
	    << setw(10) << "NA" << "\n";

      MHOUT << setw(4) << locus[l]->chr << " " 
	    << setw(par::pp_maxsnp) << locus[l]->name << " "
	    << setw(4) << locus[l]->allele1 << " "
	    << setw(4) << locus[l]->allele2 << " "
	    << setw(8) << "NA" << " "
	    << setw(8) << "NA" << " "
	    << setw(8) << "NA" << " "
	    << setw(8) << "NA" << " "
	    << setw(6) << "ASSOC" << " "
	    << setw(10) << X_assoc << " "
	    << setw(4) << 1 << " "
	    << setw(10) << chiprobP(X_assoc,1) << " "
	    << setw(10) << "NA" << "\n";

      MHOUT << setw(4) << locus[l]->chr << " " 
	    << setw(par::pp_maxsnp) << locus[l]->name << " "
	    << setw(4) << locus[l]->allele1 << " "
	    << setw(4) << locus[l]->allele2 << " "
	    << setw(8) << "NA" << " "
	    << setw(8) << "NA" << " "
	    << setw(8) << "NA" << " "
	    << setw(8) << "NA" << " "
	    << setw(6) << "HOMOG" << " "
	    << setw(10) << X_homog << " "
	    << setw(4) << nk-1 << " "
	    << setw(10) << chiprobP(X_homog,nk-1) << " "
	    << setw(10) << "NA" << "\n";

      for (int k=0; k<nk; k++)
	{

	  if ( n_11[k] + n_12[k] <= 1.0001 || 
	       n_21[k] + n_22[k] <= 1.0001 ) 
	    {

	      MHOUT << setw(4) << locus[l]->chr << " " 
		    << setw(par::pp_maxsnp) << locus[l]->name << " "
		    << setw(4) << locus[l]->allele1 << " "
		    << setw(4) << locus[l]->allele2 << " "
		    << setw(8) << "NA" << " "
		    << setw(8) << "NA" << " "
		    << setw(8) << n_11[k] + n_12[k] - 1 << " "
		    << setw(8) << n_21[k] + n_22[k] - 1 << " "
		    << setw(6) << kname[k] << " "
		    << setw(10) << "NA" << " "
		    << setw(4) << "NA" << " "
		    << setw(10) << "NA" << " "
		    << setw(10) << "NA" << "\n";
	    }
	  else
	    {
	      
	      MHOUT << setw(4) << locus[l]->chr << " " 
		    << setw(par::pp_maxsnp) << locus[l]->name << " "
		    << setw(4) << locus[l]->allele1 << " "
		    << setw(4) << locus[l]->allele2 << " "
		    << setw(8) << n_11[k]/double(n_11[k]+n_12[k]) << " "
		    << setw(8) << n_21[k]/double(n_21[k]+n_22[k]) << " "
		    << setw(8) << n_11[k] + n_12[k] - 1 << " "
		    << setw(8) << n_21[k] + n_22[k] - 1 << " "
		    << setw(6) << kname[k] << " "
		    << setw(10) << X_indiv[k] << " "
		    << setw(4) << 1 << " "
		    << setw(10) << chiprobP(X_indiv[k],1) <<  " ";
	      double odr = ( n_11[k] * n_22[k] ) / ( n_12[k] * n_21[k] );
	      if ( realnum(odr) )
		MHOUT << setw(10) 
		      << odr << "\n";
	      else
		MHOUT << setw(10) 
		      << "NA" << "\n";
	      
	    }
	}

      
      // Next locus
      s++;
      l++;

    } 

  MHOUT.close();


}




vector<double> Plink::calcMantelHaenszel_ORD(vector<int> & X,
					     vector<int> & Y,
					     vector<int> & Z)
{
  
  // X is SNP coding
  // Y is phenotype (assumed to be ordinal, integers)
  // Z is cluster
  
  if (X.size() != Y.size() || Y.size() != Z.size() || X.size() != Z.size() )
    error("Internal problem:\n  problem in calcMantelHaenszel_ORD()...uneven input columns");
  
  // Determine unique elements
  int nx=0, ny=0, nz=0;
  map<int,int> mx;
  map<int,int> my;
  map<int,int> mz;
  
  for (unsigned int i=0; i<X.size(); i++)
    {
      if (mx.find(X[i]) == mx.end())
	mx.insert(make_pair(X[i],nx++));
      
      if (my.find(Y[i]) == my.end())
	my.insert(make_pair(Y[i],ny++));
      
      if (mz.find(Z[i]) == mz.end())
	mz.insert(make_pair(Z[i],nz++));
    }
  
  
  // Generic function to calculate generalized ordinal IxJxK CMH
  // Assumes no missing data
  
  vector< vector<double> > N(nz);             // observed counts
  vector< vector<double> > U(nz);             // expected
  vector< vector< vector<double> > > V(nz);   // variance matrix
  
  vector<vector<double> > Tx(nz); // marginal totals
  vector<vector<double> > Ty(nz); // ..
  vector<double> T(nz);           // totals (per K)
  
  for (int k=0; k<nz; k++)
    {
      Tx[k].resize(nx);
      Ty[k].resize(ny);
      
      N[k].resize((nx-1)*(ny-1));
      U[k].resize((nx-1)*(ny-1));
      V[k].resize((nx-1)*(ny-1));
      for (int k2=0; k2<(nx-1)*(ny-1); k2++)
	{
	  N[k][k2] = U[k][k2] = 0;
	  V[k][k2].resize((nx-1)*(ny-1));
	  for (int k3=0; k3<(nx-1)*(ny-1); k3++)
	    V[k][k2][k3] = 0;
	}
    }
  
  // Create counts
  // Consider each observation
  for (int i=0; i<X.size(); i++)
    {
      int vx = mx.find(X[i])->second;
      int vy = my.find(Y[i])->second;
      int vz = mz.find(Z[i])->second;
      
      // exclude nx + ny (upper limits)
      if (vx<nx-1 && vy<ny-1)
	N[vz][ vx + vy*(nx-1) ]++;
      
      Tx[vz][vx]++;
      Ty[vz][vy]++;
      T[vz]++;
    }
  
  
  // Determine valid clusters (at least 2 people)
  vector<bool> validK(nk,false);
  for (int k=0; k<nk; k++)
    if (T[k]>=2) validK[k]=true;
  

  // Calculate expecteds
  for (int k=0; k<nz; k++)
    {
      if (validK[k])
	{
	  for (int ix=0; ix<nx-1; ix++)
	    for (int iy=0; iy<ny-1; iy++)
	      {
		U[k][ix+iy*(nx-1)] = ( Tx[k][ix] * Ty[k][iy] ) / T[k];
	      }
	}
    }


  

  ////////////////////
  // Ordinal CMH test
  
  // ordered scores u_i and v_i
  
// T_k = sum_i sum_j u_i v_j n_ijk


// E(T_k) = [ ( \sum_i u_i n_{i+k} ) ( \sum_j v_j n_{+jk} ) ] / n_{++k}

// var(T_k) = ( 1 / ( n_{++k} - 1 )  
//       * [ \sum_i u_i^2  n_{i+k} - ( ( \sum_i u_i n_{i+k} )^2 / n_{++k} ) ]
//       * [ \sum_j v_j^2  n_{+jk} - ( ( \sum_j v_j n_{+jk} )^2 / n_{++k} ) ]

// statistic [ T_k - E(T_k)]/[var(T_k)]^{1/2} - correl between X and Y
// in stratum k, multiplied by sqrt( n_{++k} - 1 ) 

// To summarize across K strata: 

// M^2 = ( \sum_k  [ T_k - E(T_k) ] )^2  / sum_k var(T_k) 

// which has chi-sq 1 df null distrib

  
  
}




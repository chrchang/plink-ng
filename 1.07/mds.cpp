


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
#include <fstream>
#include <cmath>

#include "plink.h"
#include "options.h"
#include "helper.h"
#include "stats.h"

#ifdef WITH_LAPACK
#include "lapackf.h"
#endif


void Plink::generateMDS()

{
  
  // Take this solution, (i)
  // 1) Average clusters to generate points
  // 2) Perform multidimensional scaling
  // 3) Dump information into Haploview-friendly file for visualisation
  
  string f = par::output_file_name + ".mds";
  printLOG("Writing MDS solution to [ " + f + " ] \n");
  
  if (par::mds_by_individual)
    printLOG("MDS plot of individuals (not clusters)\n");
  else
    printLOG("MDS plot of clusters (not individuals)\n");
  
  vector< vector<int> > cl;

  // Re-Populate the cl cluster information list, if need be
  if ( ! par::mds_by_individual )
    {
      set<int> clnum;
      for (int i=0; i<n; i++)
	{
	  if ( sample[i]->sol >= 0 )
	    clnum.insert( sample[i]->sol );
	}
      cl.resize( clnum.size() );
      for (int i=0; i<n; i++)      
	{
	  if ( sample[i]->sol >= 0 )
	    cl[ sample[i]->sol ].push_back(i);
	}
    }
 
  int nc;
  
  if (par::mds_by_individual)
    nc = n;
  else
    nc = cl.size();


  // Now we have built the between-cluster distance matrix (which will
  // typically be smaller than the between-individual matrix, we
  // should be able to apply visualisation (note: for samples of under
  // 5000 individuals, should be okay to apply standard per-individual
  // clustering
  
  // B = - 1/2 Z D^2 Z 
  // 
  // where Z = I - 1/n U 
  //
  // I identity matrix (n x n)
  // U is unit matix (n x n)
  
  // A double-centered matrix B
  // b_ij = -1/2 [ d^2_ij - d^2_.j - d^2_i. + d^2_.. ] 


#ifdef WITH_LAPACK

  // Full, symm matrix (1D format)
  vector_t D(nc*nc,0);
  
  for (int c1 = 0 ; c1 < nc ; c1++)
    for (int c2 = c1 ; c2 < nc ; c2++)
      {
	if (c1==c2) D[ c1 + c2*nc ]=0;
	else
	  {
	    if (par::mds_by_individual)
	      {
		if (c1>c2)
		  D[c1 + c2*nc] = D[ c2 + c1*nc ] = (1-mdist[c1][c2]) * (1-mdist[c1][c2]);
		else
		  D[c1 + c2*nc] = D[ c2 + c1*nc ] = (1-mdist[c2][c1]) * (1-mdist[c2][c1]);
	      }
	    else
	      {
		// Average over all pairs between cluster
		
		double avg = 0;
		for (int i1=0; i1<cl[c1].size(); i1++)
		  for (int i2=0; i2<cl[c2].size(); i2++)
		    {
		      if ( cl[c1][i1] > cl[c2][i2] ) 
			avg += 1-mdist[cl[c1][i1]][cl[c2][i2]];
		      else
			avg += 1-mdist[cl[c2][i2]][cl[c1][i1]];
		    }
		avg /= cl[c1].size() * cl[c2].size();
		
		// Symmetric matrix of squared distances
		D[c1 + c2*nc] = D[c2 + c1*nc] = avg * avg;
	      }
	    
	  }
      }
    
  double mean = 0;
  vector_t M(nc,0);
  
  for (int c1 = 0 ; c1 < nc ; c1++)
    {
      for (int c2 = 0 ; c2 < nc ; c2++)
	{
	  M[c1] += D[c1 + c2*nc];
	}
      M[c1] /= (double)nc;
      mean += M[c1];
    }
  mean /= (double)nc;
  
  // For each element for D, double center
  
  for (int c1 = 0 ; c1 < nc ; c1++)
    for (int c2 = c1 ; c2 < nc ; c2++)
      D[c1 + c2*nc] = D[c2 + c1*nc] = - 0.5 * ( D[c1 + c2*nc] - M[c1] - M[c2] + mean );

  
  // Calculate only required eigen-vectors
  vector_t eigenvalue(nc);
  matrix_t eigenvector;
  sizeMatrix(eigenvector,nc,nc);
  //svd_lapack(n,D,eigenvalue,eigenvector);
  eigen_lapack(n,D,eigenvalue,eigenvector);
  
//   cout << "EVAL = \n";
//   display(eigenvalue);
//   cout << "EVEC = \n";
//   display(eigenvector);

#else

  
  // Not using LAPACK

  matrix_t D;
  sizeMatrix(D,nc,nc);
  
  for (int c1 = 0 ; c1 < nc ; c1++)
    for (int c2 = c1 ; c2 < nc ; c2++)
      {
	if (c1==c2) D[c1][c2]=0;
	else
	  {
	    if (par::mds_by_individual)
	      {
		if (c1>c2)
		  D[c1][c2] = D[c2][c1] = (1-mdist[c1][c2]) * (1-mdist[c1][c2]);
		else
		  D[c1][c2] = D[c2][c1] = (1-mdist[c2][c1]) * (1-mdist[c2][c1]);
	      }
	    else
	      {
		// Average over all pairs between cluster
		
		double avg = 0;
		for (int i1=0; i1<cl[c1].size(); i1++)
		  for (int i2=0; i2<cl[c2].size(); i2++)
		    {
		      if ( cl[c1][i1] > cl[c2][i2] ) 
			avg += 1-mdist[cl[c1][i1]][cl[c2][i2]];
		      else
			avg += 1-mdist[cl[c2][i2]][cl[c1][i1]];
		    }
		avg /= cl[c1].size() * cl[c2].size();
		
		// Symmetric matrix of squared distances
		D[c1][c2] = D[c2][c1] = avg * avg;
	      }
	    
	  }
      }
  
  
  double mean = 0;
  vector_t M(nc,0);
  
  for (int c1 = 0 ; c1 < nc ; c1++)
    {
      for (int c2 = 0 ; c2 < nc ; c2++)
	{
	  M[c1] += D[c1][c2];
	}
      M[c1] /= (double)nc;
      mean += M[c1];
    }
  mean /= (double)nc;
  
  // For each element for D, double center
  
  for (int c1 = 0 ; c1 < nc ; c1++)
    for (int c2 = c1 ; c2 < nc ; c2++)
      D[c1][c2] = D[c2][c1] = - 0.5 * ( D[c1][c2] - M[c1] - M[c2] + mean );  

  vector_t eigenvalue(nc);
  matrix_t eigenvector;
  sizeMatrix(eigenvector,nc,nc);


//   cout << "*---------\n";
//   for (int i=0; i<nc; i++)
//     {
//       for (int j=0; j<nc; j++)
// 	cout << D[i][j] << "  ";
//       cout << "\n";
//     }
//   cout << "*---------\n";

  bool flag = svd(D,eigenvalue,eigenvector); 

//   cout << "EVAL = \n";
//   display(eigenvalue);
//   cout << "EVEC = \n";
//   display(eigenvector);

#endif
  



  /////////////////////////////////////////////////////
  // Done all SVD calculation, return to normal code


  // Take the e largest eignevectors
  map<double,int> emap;
  for (int i=0; i<nc; i++)
    emap.insert(make_pair( eigenvalue[i] , i ) );
  
  map<double,int>::reverse_iterator e = emap.rbegin();
  int inc = par::cluster_mds_dim;

  vector<int> elist;
  while ( e != emap.rend() && inc > 0 ) 
    {
      elist.push_back(e->second);
      inc--;
      e++;
    }

  if (par::cluster_mds_dim < 1) 
    par::cluster_mds_dim = 1;
  if (par::cluster_mds_dim > nc) 
    par::cluster_mds_dim = nc;
  
  if ( elist.size() != par::cluster_mds_dim )
    {
      error("Internal problem extracting MDS solution\n");
      elist.resize(par::cluster_mds_dim);
    }
  
  
  // Sqrt(D)

  for (int i=0; i<nc; i++) 
    eigenvalue[i] = eigenvalue[i] >= 0 ? sqrt(eigenvalue[i]) : 0 ;

  // Make solution
  // EVEC * sqrt(EVAL)  but filter on rows that are in solution
  // with EVAL as diagonal matrix

  matrix_t mds;
  sizeMatrix(mds,nc,par::cluster_mds_dim);

  for (int c1=0; c1<nc; c1++)
    for (int c2=0; c2<par::cluster_mds_dim; c2++)
      // ERROR: *** in 1.02 and below was ***   for (int c3=0; c3<nc; c3++)
      for (int c3=0; c3<par::cluster_mds_dim; c3++)
	{
	  int i2 = elist[c2];
	  int i3 = elist[c3];
	  if ( i3 == i2 )
 	    mds[c1][c2] += eigenvector[c1][i3] * eigenvalue[i2];
	}

  
  // Display solution
  ofstream MDS(f.c_str(),ios::out);
  MDS.precision(6);
  
  MDS << setw(par::pp_maxfid) << "FID" << " "
      << setw(par::pp_maxiid) << "IID" << " "
      << setw(6) << "SOL" << " ";
  
  for (int c=0; c<par::cluster_mds_dim; c++)
    MDS << setw(12) << "C"+int2str(c+1) << " ";
  MDS << "\n";
  
  for (int i=0; i<n; i++)
    {
      MDS << setw(par::pp_maxfid) << sample[i]->fid << " "
	  << setw(par::pp_maxiid) << sample[i]->iid << " "
	  << setw(6) << sample[i]->sol << " ";
      
      for (int c=0; c<par::cluster_mds_dim; c++)
	{
	  if ( par::mds_by_individual )
	    MDS << setw(12) << mds[i][c] << " ";
	  else
	    {
	      if ( sample[i]->sol >= 0 )
		MDS << setw(12) << mds[ sample[i]->sol ][c] << " ";
	      else
		MDS << setw(12) << "NA" << " ";
	    }
	}

      MDS << "\n";           
    }
  MDS.close();  

}

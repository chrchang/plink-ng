

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2008 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////



#ifndef __STATS_H__
#define __STATS_H__

#include <string>
#include <vector>
#include <cstdio>

#include "plink.h"

using namespace std;

void sizeMatrix(matrix_t &, int,int);
void sizeTable(table_t & , int, int);
void multMatrix(matrix_t & a,
		matrix_t & b,
		matrix_t & c);
matrix_t vec2diag(vector_t &);

class Eigen
{
 public:
  void set(int n)
    {
      d.resize(n,0);
      sizeMatrix(z,n,n);
    }
  
  vector_t d; // eigenvalues
  matrix_t z; // eigenvectors
};


bool realnum(double);

long double factorial(int);
double normdist(double);
double ltqnorm(double);
double chi2x2(double,double,double,double);
double chi2x2(table_t);
double chi2x2(matrix_t);
double chiTable(table_t);
double chiprobP(double, double);
double symTable(table_t);
double inverse_chiprob(double, double); 
double gammp(double a, double x);
void gser(double *gamser, double a, double x, double *gln);
void gcf(double *gammcf, double a, double x, double *gln);
double gammln(double xx);

double rnorm();

void lubksb(vector<vector<double> > &a, vector<int> &indx, vector<double> &b);
void ludcmp(vector<vector<double> > &a, vector<int> &indx, double &d);
vector< vector<double> > inverse(vector< vector<double> > & m );

vector<double> eigenvalues(vector<vector<double> > & a);
void tred2(vector<vector<double> >&,vector<double> &,vector<double> &);
void tqli(vector<double> &d, vector<double>&e, vector<vector<double> > &z);

Eigen eigenvectors(vector<vector<double> > & a);
void EV_tred2(vector<vector<double> >&,vector<double> &,vector<double> &);
void EV_tqli(vector<double> &d, vector<double>&e, vector<vector<double> > &z);

vector< vector<double> > svd_inverse(vector< vector<double> > & , bool & );
bool svd(matrix_t &,vector_t &, matrix_t &);
bool svdcmp(vector<vector<double> > &, 
	    vector<double> &, 
	    vector<vector<double> > &);
void svbksb(vector<vector<double> > &u,
	    vector<double> &w,
	    vector<vector<double> > &v,
	    vector<double> &b, 
	    vector<double> &x);
vector<vector<double> > msqrt(vector<vector<double> > & u);

double qromb(double func(const double), double a, double b);
void polint(vector_t &xa, vector_t &ya, const double x, double &y, double &dy);
double trapzd(double func(const double), const double a, const double b, const int n);

void svdvar(vector<vector<double> > & v,
	    vector<double> & w,
	    vector<vector<double> > & cvm);

int pca(matrix_t & x, boolmatrix_t & mask, vector_t & p, matrix_t & s,matrix_t & v, bool);

double pythag(const double a, const double b);

double betacf(const double a, const double b, const double x);
double betai(const double a, const double b, const double x);
double pF(const double F, const int df1, const int df2);
double pT(const double T, const double df);

#endif

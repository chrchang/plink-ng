

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2008 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#include <iostream>
#include <sstream>
#include <iomanip>

#include "plink.h"
#include "helper.h"
#include "options.h"
#include "stats.h"

using namespace std;

#ifdef WITH_R_PLUGINS
#define MAIN         // we are the main program, we need to define this
#define SOCK_ERRORS  // we will use verbose socket errors
#include "sisocks.h"
#include "Rconnection.h"
#endif



void Plink::Rfunc() 
{

  
#ifdef WITH_R_PLUGINS

  bool write_script = par::run_R_write_script;


  // Ensure SNP-major mode

  if ( ! par::SNP_major ) 
    Ind2SNP();

  // Are thre individuals / SNPs worth testing

  if ( n == 0 ) 
    error("No individuals left for analysis");
  else if ( nl_all == 0 ) 
    error("No SNPs left for analysis");

#ifdef SKIP
  printLOG("R-extensions not implemented on this system...\n");
  return;
#else

  printLOG("R-extension call for script [ " + par::R_script + " ]\n");
  
  checkFileExists( par::R_script ); 
  
  ifstream RIN(par::R_script.c_str(), ios::in);
  ofstream ROUT;
  ofstream RSCRIPT;

  Rconnection *rc;

  if ( ! write_script ) 
    {
      
      ROUT.open((par::output_file_name+".auto.R").c_str(), ios::out);

      printLOG("Writing results of R-extension to [ " 
	       + par::output_file_name+".auto.R ]\n");      
      
      rc = new Rconnection("127.0.0.1", par::R_port);

      
      int i=rc->connect();
      
      if (i) {
	char msg[128];
	sockerrorchecks(msg, 128, -1);
	printf("unable to connect (result=%d, socket:%s).\n", i, msg); 
      }
 
      // Minimal output    
      rc->eval("options(echo=F)");
      
      
   }
  else
    {
      printLOG("Writing debug-mode R-extension to [ " 
	       + par::output_file_name+".debug.R ]\n");      

      RSCRIPT.open((par::output_file_name+".debug.R").c_str(), ios::out);
    }
    



  /////////////////////////////////////////////
  // Read R script that defines function Rplink 

  string rcommand_data, user_function, line;
  
  while(getline(RIN, line))
    user_function += line + "\n";




  ////////////////////////////////////
  // create R friendly data structures

  // Remove individuals with missing phenotypes

  removeMissingPhenotypes(*this);


  // We are passing 'n' individuals 
  // and nl_all SNPs.  By default, we will pass all individuals, 
  // but run_R_nsnps-sized batches of SNPs only 
  
  // Phenotypes
  vector_t pvec(n);
  double * p = &(pvec[0]);
  for (int i=0;i<n;i++)
    p[i] = sample[i]->phenotype;
  
  
  
  // Cluster information, per person
  
  vector<int> svec(n);
  int * s = &(svec[0]);
  for ( int i = 0; i < n; i++ )
    s[i] = sample[i]->sol;


  // Covariates
  int x = -1;
  int size = n * par::clist_number;
  vector_t cvec(size);  
  double * c = &(cvec[0]);

  // if there exists a covariate matrix

  if( par::clist_number > 0 ){
    for( int i = 0; i < n; i++ ){
      for( int j = 0; j < par::clist_number; j++ ){
	c[++x] = sample[i]->clist[j];
      }
    }
  }

  // Assign space
  
  Rinteger * rN;
  Rdouble * rP;
  Rinteger * rS;
  Rdouble * rCov;
  
  // Assign variables in R

  if ( ! write_script ) 
    {
      
      rS = new Rinteger(s, n);
      rP = new Rdouble(p, n);
      rN = new Rinteger(&n, 1);
      
      rc->assign("n", rN );      
      rc->assign("PHENO", rP );
      rc->assign("CLUSTER", rS);
      
      rc->eval("CLUSTER[CLUSTER==-1] <- NA");

      // If there exists a covariate matrix

      if( par::clist_number > 0 ){
	rCov =  new Rdouble(c, size);
	rc->assign("c", rCov);
	rc->eval("COVAR<-matrix(c,nrow=n,byrow=T)");
      }
      else{
	rc->eval("COVAR<-NA");
      }
      
    }
  else
    {
      
      // Write commands to file
      
      RSCRIPT << "n <- " << n << "\n";
      
      RSCRIPT << "PHENO <- c( ";
      for (int i=0; i<n-1; i++)
	RSCRIPT << p[i] << ", ";
      RSCRIPT << p[n-1] << " ) \n";

      if ( par::clist_number > 0 && n > 0)
	{
	  RSCRIPT << "c <- c( ";
	  for (int i=0; i<n*par::clist_number-1; i++)
	    RSCRIPT << c[i] << ", ";
	  RSCRIPT << c[n*par::clist_number-1] << " ) \n";
	  RSCRIPT << "COVAR <- matrix( c , nrow = n , byrow=T)\n";	  
	}
      else
	{
	  RSCRIPT << "COVAR <- matrix( NA , nrow = n , ncol = 0 , byrow = T)\n";	  
	}

      RSCRIPT << "CLUSTER <- c( ";
      for (int i=0; i<n-1; i++)
	RSCRIPT << s[i] << ", ";
      RSCRIPT << s[n-1] << " ) \n";
      
    }
   

  // Now loop over batches of genotypes, in batches

  for ( int l=0 ; l < nl_all ; l += par::run_R_nsnps )
    {
      
      int nstart = l;
      int nstop = ( nl_all-1 ) < l + par::run_R_nsnps - 1 ?
	nl_all-1 : l + par::run_R_nsnps - 1 ;
      int nloc = nstop - nstart + 1;
      
      if ( ! ( par::silent || write_script ) )
	{
	  cout << "Considering SNPs from " 
	       << nstart 
	       << " to " << nstop 
	       << "             \r";
	  cout.flush();
	}

      // Genotypes ( # of minor allele)
      
      vector<int> gvec(n*nloc);
      int * g = &(gvec[0]);
      int x = 0;
          
      for ( int i = 0; i < n; i++ )
	for( int j = nstart; j <= nstop; j++ ){
	  
	  bool one = SNP[j]->one[i];
	  bool two = SNP[j]->two[i];
	  
	  if( (!one) && !two)
	    g[x] = 2;
	  if( (!one) && two )
	    g[x] = 1;
	  if( one && (!two) )
	    g[x] = -1;
	  if( one && two)
	    g[x] = 0;
	  
	  x++;
	}
      
      Rinteger * rL = new Rinteger(&nloc,1);
      Rinteger * rG = new Rinteger(g, n*nloc);

      // Assign variables in R
            
      if ( ! write_script ) 
	{
	  rc->assign("l", rL );
	  rc->assign("g", rG );
     
	  rc->eval("GENO<-matrix(g,nrow=n,byrow=T)");
	  rc->eval("GENO[GENO==-1] <- NA");
	}
      else
	{
	  // Write commands to file
	  
	  RSCRIPT << "l <- " << nloc << "\n";
	  
	  RSCRIPT << "g <- c( ";
	  for (int i=0; i<n*nloc-1; i++)
	    RSCRIPT << g[i] << ", ";
	  RSCRIPT << g[n*nloc-1] << " ) \n";
	  RSCRIPT << "GENO <- matrix( g , nrow = n ,byrow=T)\n";
	  RSCRIPT << "GENO[GENO == -1 ] <- NA \n";
	}
      
      // free space
      delete rL;
      delete rG;  
      
      if ( write_script ) 
	{
	  RSCRIPT << "\n\n" << user_function << "\n";
	  continue;
	}

      //////////////////////////////////////////////////
      // Run R script which will create function Rplink
      
      rc->eval( user_function.c_str() );

      ///////////////////////////////////////////////////////////
      // And call the user's function, saving vector of results
      
      Rdouble *data = (Rdouble*) rc->eval("Rplink(PHENO,GENO,CLUSTER,COVAR)");

      ///////////////////////////////////////////////////////////
      // If everything went okay, we can get the results
      
      if (data) 
	{ 
	  // Store the results in a vector of doubles
	  
	  double * d = data->doubleArray();
	  
	  // expect format
	  // N { N items per SNP } 
	  
	  int i = 0; 
	  int ct = data->length();
		  
	  for (int l=nstart; l < nstart + nloc; l++)
	    {
	      
	      ROUT << setw(4) << locus[l]->chr << " "
		   << setw(par::pp_maxsnp) << locus[l]->name << " "
		   << setw(10) << locus[l]->bp << " "
		   << setw(4) << locus[l]->allele1 << " ";
	      
	      int c = (int)d[i++];
	      
	      for (int j=0;j<c;j++)
		{
		  if ( realnum( d[i] ) )
		    ROUT << d[i++] << "\t";
		  else
		    {
		      ROUT << "NA" << "\t";
		      ++i;
		    }
		}
	      ROUT << "\n";
	    }

	  // Dispose of the object
	  delete data;
	}
      else
	{

	  // populate results with missing code

	  for (int l=nstart; l < nstart + nloc; l++)
	    ROUT << setw(4) << locus[l]->chr << " "
		 << setw(par::pp_maxsnp) << locus[l]->name << " "
		 << setw(10) << locus[l]->bp << " "
		 << setw(4) << locus[l]->allele1 << " "
		 << "NA" << "\n";	  
	}
      
    }
  

  if ( ! write_script ) 
    {
      delete rN;
      delete rP;
      delete rS;
      if( par::clist_number > 0 )
	delete rCov;
    }
  
  // If in DEBUG mode, now close this file

  if ( write_script ) 
    {
      RSCRIPT.close();
      shutdown();
    }
  

  if ( ! par::silent ) 
    cout << "\n";
  
  
  // Dispose the connection object, which implicitly closes the
  // connection
  
  delete rc;


  ROUT.close();
   
  
#endif
#endif
 

  return;
  
}

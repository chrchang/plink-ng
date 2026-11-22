

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2008 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#include "gvar.h"
#include "helper.h"
#include "options.h"
#include "plink.h"
#include "model.h"
#include "logistic.h"
#include "linear.h"
#include "stats.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <errno.h>


// A couple of helper functions
Model * analyseModel(Plink *, Variant *, int, bool, bool);
double compareModels(Model *, Model *) ;


class fullGenotype{
public:
  fullGenotype(int2 a, int2 b)
  {
    a1 = a;
    a2 = b;
    if ( a1 < a2 ) 
      {
 	int2 t = a1;
 	a1 = a2;
 	a2 = t;
      }
  }
  int2 a1;
  int2 a2;

  bool operator< (const fullGenotype & b) const
  {
    return (a1 < b.a1 || (a1 == b.a1 && a2 < b.a2) );
  }
  
};

void Plink::readGenericVariantData()
{
  
  ///////////////////////////////////////////////
  // Load data (if par::gvar is true) -- otherwise, 
  // we are just
    
  string filename = par::gvarfile;

  if ( par::gvar )
    checkFileExists( filename );


  // If we are not reading in any generic variants, 
  // but we have elected to write a file, then 
  // we must set this flag so that the data are 
  // carried over

  if ( par::gvar_write && ! par::gvar ) 
    par::gvar_include_all_variants = true;
   
  ///////////////////////////////////////////////
  // Input mode
  //
  // 1) Either we start from scratch, read in gMAP file 
  //    and then data (par::load_gvar is true)

  // 2) Or, we have first read in the standard file, and we 
  //    just place these on top, copying basic SNPs over first
  //    (par::load_gvar is false)

  bool preexisting = locus.size() > 0;


  // We need either basic SNPs or generic variants here

  if ( ! ( preexisting || par::gvar ) )
    return;


  if ( preexisting )
    {
      nl_all = locus.size();
      n = sample.size();
      // prettyPrintLengths();
    }


  ///////////////////////////////////////////////
  // .map file
  
  vector<bool> include;
  vector<int> include_pos(0);

  if ( par::gvar && !preexisting )
    {

      checkFileExists( par::gmapfile );

      ngvar = 0;
      readMapFile(par::gmapfile,
		  include,
		  include_pos,
		  ngvar);
      
      gvar.clear(); 
      for (int l=0; l<locus.size(); l++)
	{
	  gvar.push_back( (Variant *)(locus[l]) );
	  gvar[l]->acode.clear();
	}
      locus.clear();

    }
  

  ///////////////////////////////////////////////
  // Otherwise figure out which markers to keep
  // i.e. assume a small subset of all

  set<string> observedGVARs;
  map<int,int> copyover;

  if ( preexisting )
    {

      // Either take all SNPs in memory, or 
      // select out just those already specified in 
      // the GVAR file

      if ( par::gvar_include_all_variants )
	{
	  // Copy all SNPs over
	  for (int l=0; l<nl_all; l++)
	    observedGVARs.insert( locus[l]->name );
	  ngvar = observedGVARs.size();

	}
      else if ( par::gvar )
	{
	  FILE * GV;      
	  GV = fopen64( filename.c_str(),"r");
	  if ( GV == NULL )
	    error("Problem opening GVAR file, errno = "+int2str(errno));
	  while ( ! feof(GV) )
	    {
	      string dummy;
	      string gvarname;
	      int f=0;
	      if ( readString( GV , dummy ) ) ++f;
	      if ( dummy == "" ) 
		continue;
	      if ( readString( GV , dummy ) ) ++f;
	      if ( readString( GV , gvarname ) ) 
		observedGVARs.insert( gvarname );
	      while (fgetc(GV) != '\n' && !feof(GV)) {}
	    }
	  
	  fclose(GV);
	  
	  ngvar = observedGVARs.size();
	  printLOG("Found "+int2str( ngvar )+" variant in [ " + filename + " ]\n");

	}


      

      //////////////////////////////////
      // Copy MAP into appropriate space

      gvar.clear(); 
      map<string,int> im;
      for (int l=0; l<locus.size(); l++)
	im.insert(make_pair( locus[l]->name, l));

      set<string>::iterator is = observedGVARs.begin();
      while ( is != observedGVARs.end() )
	{
	  map<string,int>::iterator imi = im.find( *is );
	  if ( imi != im.end() )
	    {
	      Variant * gloc = new Variant;
	      Locus * loc = locus[imi->second];
	      gloc->name = loc->name;
	      gloc->chr = loc->chr;
	      gloc->bp = loc->bp;
	      gloc->pos = loc->pos;
	      gloc->acode.clear();
	      gloc->alleles.clear();
	      gvar.push_back( gloc );

	      int l = gvar.size() - 1;
	      copyover.insert(make_pair(l,imi->second));
	    }
	  else
	    error("Generic variant " 
		  + *is 
		  + " found that was not in original SNP file\n");

	  ++is;
	}
    }



  ///////////////////////////////////////////////
  // .fam
  
  if ( par::gvar && ! preexisting )
    {
      checkFileExists( par::gfamfile );
      readFamFile(par::gfamfile);
      n = sample.size();

      // Read an alternate phenotype file? 
      if (par::pheno_file) 
	readPhenoFile();

      if ( par::bt )
	affCoding(*this);
    }



  ///////////////////////////////////////////////
  // Allocate space

  for (int i=0; i<sample.size(); i++)
    {
      sample[i]->gvar.resize( ngvar );
      for (int g=0; g<ngvar; g++)
	sample[i]->gvar[g] = new GVariant;
    }



  ///////////////////////////////////////////////
  // Copy over any existing SNP data

  if ( preexisting )
    {
      for (int g=0; g<ngvar; g++)
	{
	  Variant * gloc = gvar[g];
	  map<int,int>::iterator li = copyover.find(g);
	  int l = li->second;

	  Locus * loc = locus[ li->second ];
	  
	  int acode1 = -1;
	  int acode2 = -1;

	  if ( loc->allele1 != "" && loc->allele1 != par::missing_genotype )
	    {
	      if ( gloc->acode.find( loc->allele1 ) == gloc->acode.end() )
		{
		  gloc->acode.insert( make_pair( loc->allele1, gloc->nallele ));
		  gloc->alleles.push_back( loc->allele1 );
		  acode1 = gloc->nallele;
		  ++(gloc->nallele);
		}
	    }

	  if ( loc->allele2 != "" && loc->allele2 != par::missing_genotype )
	    {
	      if ( gloc->acode.find( loc->allele2 ) == gloc->acode.end() )
		{
		  gloc->acode.insert( make_pair( loc->allele2, gloc->nallele ));
		  gloc->alleles.push_back( loc->allele2 );
		  acode2 = gloc->nallele;
		  ++(gloc->nallele);
		}
	    }


	  // Now copy over actual genotypes for this particular SNP
	  // allele

	  for (int i=0; i<sample.size(); i++)
	    {

	      GVariant * gv = sample[i]->gvar[g];

	      gv->allele1 = gv->allele2 = -1;

	      bool s1 = par::SNP_major ? SNP[l]->one[i] : sample[i]->one[l];
	      bool s2 = par::SNP_major ? SNP[l]->two[i] : sample[i]->two[l];

	      // Missing genotype
	      if ( s1 && ! s2 )
	      {
		gv->missing = true;
		continue;
	      }
	      
	      gv->missing = false;
	      
	      if ( s1 )
		{
		  gv->allele1 = gv->allele2 = acode2;
		}
	      else
		{
		  if ( s2 )
		    {
// 		      if ( sample[i]->sex ) 
// 			cout << "het male??? " << sample[i]->fid << "\n";
		      
		      gv->allele1 = acode1;		      
		      gv->allele2 = acode2;
		    }
		  else
		    {
		      gv->allele1 = acode1;
		      gv->allele2 = acode1;
		    }
		}


	      ///////////////////////////////////////////// 	      
	      // Autosomal or haploid? Should have already 
	      // blanked out hemizygous haploid calls
	      
	      if ( par::chr_haploid[gloc->chr] || 
		   ( par::chr_sex[gloc->chr] && sample[i]->sex ) )
		{
		  gv->dosage1 = 1;
		  gv->dosage2 = 0;
		}
	      else
		{
		  gv->dosage1 = 1;
		  gv->dosage2 = 1;
		}
		
	      
	    }
	}
    }



  ///////////////////////////////////////////////
  // .gvar
  
  FILE * GV;
  
  if ( par::gvar )
    {
      GV = fopen64( filename.c_str(),"r");
      if ( GV == NULL )
	error("Problem opening GVAR file, errno = "+int2str(errno));
  
      
  // We can now read any number of individual/genotype lines, in any
  // order; we also do not assume that all genotypes are given --
  // these will be missing by default
  
  map<string,int> imap;
  if ( ! preexisting )
    {
      for (int i=0; i<include.size(); i++)
	{
	  if ( include[i] ) 
	    {
	      int k = include_pos[i];
	      imap.insert( make_pair ( locus[k]->name , k ) );
	    }
	}
    }
  else
    {
      include.clear();
      include.resize( ngvar, true );
      include_pos.clear();
      for (int i=0; i<include.size(); i++)
	{
	  include_pos.push_back(i);
	  imap.insert( make_pair ( gvar[i]->name , i ) );
	}
    }
  
  map<string,int> iperson;  
  for (int i=0; i<sample.size(); i++)
    {
      iperson.insert( make_pair 
		      ( sample[i]->fid + "_" + sample[i]->iid , 
			i ) );
    }
  
  
  // Whether or not we want to look at a locus is in the include[] vector
  // The genomic position of locus i is k=include_pos[i] -> locus[k]

  bool fatal = false;
  string fmsg = "";

  while( ! feof(GV) )
    {
      
      string fid = "";
      string iid = "";
      string gvarname = "";
      string one = "";
      string dose1 = "";
      string two = "";
      string dose2 = "";

      int f = 0;
      
      if ( readString( GV , fid ) ) f++;
      
      if ( fid == "" ) 
	continue;
      
      if ( readString( GV , iid ) ) f++;
      if ( readString( GV , gvarname ) ) f++;

      if ( readString( GV , one ) ) f++;
      if ( readString( GV , dose1 ) ) f++;

      if ( readString( GV , two ) ) f++;
      if ( readString( GV , dose2 ) ) f++;

      map<string,int>::iterator peri 
	= iperson.find( fid+"_"+iid );      
      
      Individual * person = 
	peri != iperson.end() ?
	sample[peri->second] : NULL ; 
	
      map<string,int>::iterator im = imap.find( gvarname );
      
      int k = 
	im != imap.end() ?
	im->second : -1;
      
      // Ignore this genotype?

      if ( ( ! person ) || k < 0 ) 
	continue;

      double d1, d2;

      if ( ( ! from_string<double>( d1, dose1, std::dec) ) 
	   || ( ! from_string<double>( d2, dose2, std::dec) ) )
	error("Problem in format of file [ " + filename + " ]\n\n"
	      +fid+" "+iid+" "+gvarname+" "
	      +one+" "+dose1+" "
	      +two+" "+dose2+"\n");
      
      int ip = peri->second;	
      Variant * gloc = gvar[k];
      GVariant * g = person->gvar[k];

      /////////////////////////////////////////
      // Add allele names to list, if needed

      if ( ( one == par::missing_genotype && d1 > 0 ) ||
	   ( two == par::missing_genotype && d2 > 0 ) )
	{
	  g->missing = true;
	  continue;
	}

      g->missing = false;

      map<string,int>::iterator ia = gloc->acode.find( one );


      if ( d1 > 0 )
	{
	  if ( ia == gloc->acode.end() )
	    {
	      gloc->acode.insert( make_pair( one, gloc->nallele ));
	      gloc->alleles.push_back( one );
	      g->allele1 = gloc->nallele;
	      ++(gloc->nallele);
	    }
	  else
	    g->allele1 = ia->second;
	}
      else
	g->allele1 = -1;
      

      // Store second allele

      ia = gloc->acode.find( two );

      if ( d2 > 0 )
	{
	  if ( ia == gloc->acode.end() )
	    {
	      gloc->acode.insert( make_pair( two, gloc->nallele ));
	      gloc->alleles.push_back( two );
	      g->allele2 = gloc->nallele;
	      ++(gloc->nallele);
	    }
	  else
	    g->allele2 = ia->second;
	}
      else
	g->allele2 = -1;

      // Store dosage information
      
      g->dosage1 = d1;
      g->dosage2 = d2;

      // Have we seen any non-integer dosages?

      if ( d1 - int(d1) > 1e-6 || 
	   d2 - int(d2 ) > 1e-6 ) 
	gloc->integerDosage = false;

    }

  
  fclose(GV);
}


  /////////////////////////////
  // Clear up any old SNP data
  
  if ( preexisting )
    {
      
      for (int l=0; l<locus.size(); l++)
	delete locus[l];      

      if (par::SNP_major)
	for (int l=0; l<SNP.size(); l++)
	  delete SNP[l];      
      else
	for (int i=0; i<sample.size(); i++)
	  {
	    sample[i]->one.clear();
	    sample[i]->two.clear();
	  }

      SNP.clear();
      locus.clear();
      
    }
  
}



void Plink::processGVAR()
{
  
  
  // Use the Model class to test generic variants, 
  // but not entering any main SNP effects

  par::assoc_glm_without_main_snp = true;
  
  printLOG("Processing data for " 
	   + int2str( ngvar ) 
	   + " generic variants\n");
  
  ofstream GOUT;
  ofstream GVERB;

  if ( true )
    {
      printLOG("Writing frequency & genotyping informtion to [ " 
	       + par::output_file_name 
	       + ".gvar.summary ]\n");
      GOUT.open( ( par::output_file_name 
		   + ".gvar.summary").c_str(), ios::out );
      GOUT.precision(4);

      GOUT << setw(16) << "NAME" << " "
	   << setw(12) << "FIELD" << " "
	   << setw(12) << "VALUE" << "\n";

    }
  
  if ( par::gvar_verbose_association )
    {
      GVERB.open((par::output_file_name+".assoc.gvar").c_str(),ios::out);
      printLOG("Writing verbose GVAR association results to [ " 
	       + par::output_file_name + ".assoc.gvar ]\n");
    }
  


  /////////////////////////////////////
  // Is there any phenotypic variation?

  bool phenotypicVariation = false;

  for (int i=1; i<n; i++)
    {
      if ( sample[i-1]->missing ) continue;
      if ( sample[i]->missing ) continue;
      if ( sample[i-1]->phenotype != sample[i]->phenotype )
	{
	  phenotypicVariation = true;
	  break;
	}
    }


  for (int g=0; g<ngvar; g++)
    {

      Variant * v = gvar[g];

      GOUT << setw(16) << v->name << " "
	   << setw(12) << "CHR" << " "
	   << setw(12) << v->chr << "\n";
      
      GOUT << setw(16) << v->name << " "
	   << setw(12) << "BP" << " "
	   << setw(12) << v->bp << "\n";
      
      vector_t f( v->nallele );
      
      map<int,int> cnvCount;
      int sampleTotalDose = 0;
      int sampleTotalInd = 0;
      
      map<fullGenotype,int2> typeCount;
      map<fullGenotype,int2>::iterator tit;

      for (int i=0; i<n; i++)
	{
	  Individual * person = sample[i];

	  if ( person->missing )
	    continue;

	  GVariant * gv = person->gvar[g];

	  if ( gv->missing )
	    {
	      int2 p(-1,1);
	      fullGenotype fg(p,p);
	      tit = typeCount.find(fg);
	      int2 q(0,0);
	      if ( person->aff )
		q.p1=1;
	      else
		q.p2=1;
	      if ( tit == typeCount.end() )
		typeCount.insert(make_pair(fg,q));
	      else
		{
		  if ( person->aff )
		    ++(tit->second.p1);
		  else
		    ++(tit->second.p2);
		}
	      continue;
	    }


	  /////////////////////////
	  // Calculate frequencies

	  if ( gv->allele1 >= 0 )
	    f[ gv->allele1 ] += gv->dosage1;
	  if ( gv->allele2 >= 0 )
	    f[ gv->allele2 ] += gv->dosage2;

	  int totalDose = (int)(gv->dosage1 + gv->dosage2);
	  cnvCount[totalDose]++;
	  sampleTotalDose += totalDose;

	  // Keep track of specific allele/CNV counts
	  // by case/control status
	  
	  if ( person->aff || par::qt ) 
	    {
	      int2 p(gv->allele1,(int)gv->dosage1);
	      int2 q(gv->allele2,(int)gv->dosage2);
	      fullGenotype fg(p,q);
	      tit = typeCount.find(fg);
	      if ( tit == typeCount.end() )
		typeCount.insert(make_pair(fg,int2(1,0)));
	      else
		++(tit->second.p1);
	    }
	  else
	    {
	      int2 p(gv->allele1,(int)gv->dosage1);
	      int2 q(gv->allele2,(int)gv->dosage2);
	      fullGenotype fg(p,q);
	      tit = typeCount.find(fg);
	      if ( tit == typeCount.end() )
		typeCount.insert(make_pair(fg,int2(0,1)));
	      else
		++(tit->second.p2);
	    }
	    	  
	  ++sampleTotalInd;
	
	  
	} // Look at next individual

      
      v->allelicVariation = false;
      int commonAlleles = 0;
      for ( int x=0; x< v->nallele; x++)
	{
	  f[x] /= sampleTotalDose;
	  if ( f[x] >= par::min_af )
	    ++commonAlleles;
	}

      if ( commonAlleles>1 )
	v->allelicVariation = true;
      
      v->copyNumberVariation = false;
      int commonCNV = 0;
      map<int,int>::iterator ia = cnvCount.begin();
      while ( ia != cnvCount.end() )
	{
	  if ( (double)ia->second 
	       / (double)sampleTotalInd  
	       >= par::min_af )
	    ++commonCNV;
	  ++ia;
	}
      if ( commonCNV>1 )
	v->copyNumberVariation = true;
      


      /////////////////////////////////
      // Report to summary file

      GOUT << setw(16) << v->name << " "
	   << setw(12) << "CNV" << " ";
      if ( v->copyNumberVariation )
	GOUT << setw(12) << "yes" << "\n";
      else
	GOUT << setw(12) << "no" << "\n";
      
      GOUT << setw(16) << v->name << " "
	   << setw(12) << "ALLELIC" << " ";
      if ( v->allelicVariation )
	GOUT << setw(12) << "yes" << "\n";
      else
	GOUT << setw(12) << "no" << "\n";
      
      GOUT << setw(16) << v->name << " "
	   << setw(12) << "GCOUNT" << " "
	   << setw(12) << sampleTotalInd << "\n";


      if ( v->integerDosage )
	GOUT << setw(16) << v->name << " "
	     << setw(12) << "INTEGER" << " "
	     << setw(12) << "Y" << "\n";
      else
	GOUT << setw(16) << v->name << " "
	     << setw(12) << "INTEGER" << " "
	     << setw(12) << "N" << "\n";

      for (int a = 0; a < v->nallele; a++)
	{
	  GOUT << setw(16) << v->name << " "
	       << setw(12) << v->alleles[a] << " "
	       << setw(12) << f[a] << "\n";
	}
      

      if ( v->integerDosage )
	{
	  ia = cnvCount.begin();
	  while ( ia != cnvCount.end() )
	    {
	      GOUT << setw(16) << v->name << " "
		   << setw(12) << "["+int2str(ia->first)+"]" << " "
		   << setw(12) 
		   << (double)ia->second / (double)sampleTotalInd 
		   << "\n";
	      ++ia;
	    }
	  

	  ////////////////////////////////////////////
	  // Display full Genotype counts (cases/all)
	  
	  tit = typeCount.begin();
	  while ( tit != typeCount.end() )
	    {
	      string aname = "";
	      if ( (tit->first.a1.p2 == 0 ) )
		aname += "null";
	      for (int z=0; z<tit->first.a1.p2; z++)
		{
		  if ( tit->first.a1.p1 == -1 )
		    aname += par::out_missing_genotype;
		  else
		    aname +=  v->alleles[ tit->first.a1.p1 ];
		}
	      
	      aname += "/";
	      if ( (tit->first.a2.p2 == 0 ) )
		aname += "null";
	      for (int z=0; z<tit->first.a2.p2; z++)
		{
		  if ( tit->first.a2.p1 == -1 )
		    aname += par::out_missing_genotype;
		  else
		    aname +=  v->alleles[ tit->first.a2.p1 ];
		}
	      
	      GOUT << setw(16) << v->name << " "
		   << setw(12) << aname << " ";
	      
	      if ( par::bt )
		GOUT <<  setw(12) 
		     << int2str(tit->second.p1)+":"+int2str(tit->second.p2)
		     << "\n";
	      else
		GOUT << setw(12) 
		     << tit->second.p1
		     << "\n";
	      
	      ++tit;
	    }
	  
	  
	}
      
      ///////////////////////////////
      // Test for assocaition

      if ( phenotypicVariation && 
	   ( v->allelicVariation || v->copyNumberVariation ) )
	{

	  // Place allelic variant and/or dosage information in 
	  // covariate fields
	  
	  par::clist = true;
	  
	  Model * mJoint = NULL;
	  Model * mCNV = NULL;
	  Model * mAllelic = NULL;
	  
	  
	  if ( v->allelicVariation && v->copyNumberVariation )
	    mJoint = analyseModel(this,v,g,true,true);
	  
	  if ( v->allelicVariation )
	    mAllelic = analyseModel(this,v,g,true,false);

	  if ( v->copyNumberVariation )
	    mCNV = analyseModel(this,v,g,false,true);
	  

	  ///////////////////////////////////////////////////////////////
	  // Extract results; for now we assume everything is biallelic

	  if ( mAllelic && mAllelic->isValid() ) 
	    {
	      vector_t b = mAllelic->getCoefs();
	      vector_t var = mAllelic->getVar();
	      vector_t pval = mAllelic->getPVals();

	      int term1 = par::clist_number + 1;
	      	    
//	      double chisq = mAllelic->getStatistic();
//	      double pvalue = chiprobP(chisq,1);
	      
	      GOUT << setw(16) << v->name << " "
		   << setw(12) << "B(SNP)" << " " 
		   << setw(12) << b[ term1 ] << "\n";	     
	      
	      GOUT << setw(16) << v->name << " "
		   << setw(12) << "P(SNP)" << " " 
		   << setw(12) << pval[ term1 - 1 ]  << "\n";	     
	    }
	  

	  if ( mCNV && mCNV->isValid() ) 
	    {
	      
	      vector_t b = mCNV->getCoefs();
	      vector_t var = mCNV->getVar();
	      vector_t pval = mCNV->getPVals();

	      int term1 = par::clist_number + 1;

	      // double chisq = mCNV->getStatistic();
	      // double pvalue = chiprobP(chisq,1);
	      
	      GOUT << setw(16) << v->name << " "
		   << setw(12) << "B(CNP)" << " " 
		   << setw(12) << b[ term1 ] << "\n";	     
	      
	      GOUT << setw(16) << v->name << " "
		   << setw(12) << "P(CNP)" << " " 
		   << setw(12) << pval[  term1 - 1 ] << "\n";	     
	      	      		
	    }

	  
	  if ( mJoint && mJoint->isValid() ) 
	    {
	
	      int term1 = par::clist_number + 1;
	      int term2 = par::clist_number + 2;
      
	      vector_t b = mJoint->getCoefs();	      
	      vector_t var = mJoint->getVar();
	      vector_t p = mJoint->getPVals();

	      // A 2df test (joint test of two parameters)
	      vector_t h(2,0);
	      matrix_t H; // row = number of fixes; cols = np
	      sizeMatrix(H,2,mJoint->getNP());
	      H[0][term1] = H[1][term2] = 1;
	      
	      double chisq = mJoint->isValid() ? 
		mJoint->linearHypothesis(H,h) : 0;
	      double pvalue = chiprobP(chisq,2);

 	      GOUT << setw(16) << v->name << " "
 		   << setw(12) << "B(CNP|SNP)" << " "
 		   << setw(12) << b[term1] << "\n";

 	      GOUT << setw(16) << v->name << " "
 		   << setw(12) << "P(CNP|SNP)" << " "
 		   << setw(12) << p[ term1 - 1 ] << "\n";

 	      GOUT << setw(16) << v->name << " "
 		   << setw(12) << "B(SNP|CNP)" << " "
 		   << setw(12) << b[ term2 ] << "\n";

 	      GOUT << setw(16) << v->name << " "
 		   << setw(12) << "P(SNP|CNP)" << " " 
 		   << setw(12) << p[ term2 - 1 ] << "\n";

	      GOUT << setw(16) << v->name << " "
		   << setw(12) << "P(SNP&CNP)" << " "
		   << setw(12) << pvalue << "\n";
	      
	    }



	  ///////////////////////////////////////////
	  // No valid model
	  
	  if ( ! ( mJoint || mAllelic || mCNV ) ) 
	    {
	      GOUT << setw(16) << v->name << " "
		   << setw(12) << "SNP/CNP" << " "
		   << setw(12) << "NA" << "\n";
	    }
	  
      
	  //////////////////////////////////
	  // Clean-up
	  
	  if ( mJoint )
	    delete mJoint;
	  if ( mAllelic )
	    delete mAllelic;
	  if ( mCNV )
	    delete mCNV;

  
	} // End of association testing

    } // Next generic variant
  

  if ( true ) 
    GOUT.close();
  if ( par::gvar_verbose_association )
    GVERB.close();

}





Model * analyseModel(Plink * P, Variant * v, int g, bool allelic, bool cnv)
{

  // Keep track of original numbr of covariates

  int addedTerms = 0;
  
  if ( allelic && cnv ) 
    {
      addedTerms = 2;
      P->clistname.push_back("SNP");
      P->clistname.push_back("CNP");	      
    }
  else if ( allelic ) 
    {
      addedTerms = 1;
      P->clistname.push_back("SNP");
    }
  else if ( cnv ) 
    {
      addedTerms = 1;
      P->clistname.push_back("CNP");
    }

  par::clist_number += addedTerms;

  for (int i=0; i<P->n; i++)
    {
 
      Individual * person = P->sample[i];
      GVariant * gv = person->gvar[g];
            
      // Assume biallelic for now
      
      // Y ~ m + b1.(A+B) + b2.(A-B)

      double d0=0, d1=0;

      if ( gv->missing ) 
	{
	  // This flag means that the Model 
	  // class will ignore this person

	  person->missing2 = true;
	}
      else
	{
	  person->missing2 = false;
	
	  if ( gv->allele1 == 0 )
	    d0 += gv->dosage1;
	  if ( gv->allele1 == 1 )
	    d1 += gv->dosage1;
	  
	  if ( gv->allele2 == 0 )
	    d0 += gv->dosage2;
	  if ( gv->allele2 == 1 )
	    d1 += gv->dosage2;
	}

      if ( allelic && cnv ) 
	{
	  person->clist.push_back( d0+d1 );
	  person->clist.push_back( d0-d1 );	
	}
      else if ( allelic ) 
	{
	  person->clist.push_back( d0-d1 );	
	}
      else if ( cnv ) 
	{	  
	  person->clist.push_back( d0+d1 );	
	}
	      
    }
  

  // Fit linear model, and return a pointer to it

  P->glmAssoc(false, *(P->pperm) );


  // Return covariate list to normal status

  par::clist_number -= addedTerms;
  P->clistname.resize( par::clist_number );
  for (int i=0; i<P->n; i++)
    P->sample[i]->clist.resize( par::clist_number );
  
  return P->model;
}



double compareModels(Model * alternate, Model * null)
{

  if ( par::bt )
    {      
      return chiprobP( ((LogisticModel*)null)->getLnLk() 
		       - ((LogisticModel*)alternate)->getLnLk() , 
		       alternate->getNP() - null->getNP() );
    }
  else
    {
      double F = ((LinearModel*)alternate)->calculateFTest((LinearModel*)null);
      if ( F < 0 ) F = 0;     
      return pF( F, 
		 alternate->getNP() - null->getNP(),
		 alternate->Ysize() - alternate->getNP() - 1 );           
    }
  return -1;
}


void Plink::convertGenericVariantData()
{
  error("Not yet implemented");
}

void Plink::outputGenericVariantFile()
{

  ////////////////////////////////////////////////////////////
  // Assume a fully populated generic variant dataset exists

  ////////////////////////////////
  // Write in variant-major order
  
  string f = par::output_file_name + ".gvar";
  printLOG("Writing generic variant file to [ " + f + " ]\n");

  ofstream GOUT(f.c_str(), ios::out);

  for (int g=0; g<ngvar; g++)
    {      

      Variant * v = gvar[g];
      
      for (int i=0; i<n; i++)
	{

	  Individual * person = sample[i];
	  GVariant * gv = person->gvar[g];
	  	  
	  GOUT << setw(par::pp_maxfid) << person->fid << " "
	       << setw(par::pp_maxiid) << person->iid << " "
	       << setw(par::pp_maxsnp) << v->name << " ";


	  /////////////////////
	  // First allele 

	  if ( gv->allele1 >= 0 )
	    GOUT << setw(4) << v->alleles[ gv->allele1 ] << " ";	 
	  else
	    GOUT << setw(4) << par::out_missing_genotype << " ";
	  GOUT << setw(4) << gv->dosage1 << " ";	 


	  /////////////////////
	  // Second allele 

	  if ( gv->allele2 >= 0 )
	    GOUT << setw(4) << v->alleles[ gv->allele2 ] << " ";	 
	  else
	    GOUT << setw(4) << par::out_missing_genotype << " ";
	  GOUT << setw(4) << gv->dosage2 << " ";	 
	  GOUT << "\n";

	} // Next individual
    } // Next variant

  GOUT.close();
  
}

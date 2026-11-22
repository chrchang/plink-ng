
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

#include "plink.h"
#include "helper.h"
#include "options.h"
#include "crandom.h"
#include "stats.h"
#include <cmath>

using namespace std;


//////////////////////////////////////////////////////////////////////
// A simple routine to simulate a dataset of unlinked case/control SNPs
// or QTs


class SimParameters
{
public:
  int nsnp;
  double lfreq;
  double ufreq;
  double hetOdds;
  double homOdds;
  double missing;
  string name;

  double lmarker;
  double umarker;
  double dprime;

  SimParameters()
  {
    name = "";
    nsnp = 0;
    missing = 0.00;
    lfreq = ufreq = hetOdds = homOdds = 0;
    lmarker = umarker = 0;
    dprime = 1;
  }
};


class SimParametersQT
{
public:
  int nsnp;
  double lfreq;
  double ufreq;

  double variance;
  double dom;
  double gAA, gAB, gBB;

  double missing;
  string name;

  double lmarker;
  double umarker;
  double dprime;


  SimParametersQT()
  {
    name = "";
    nsnp = 0;
    missing = 0.00;
    variance = dom = 0;
    lfreq = ufreq = 0;    
    lmarker = umarker = 0;
    dprime = 1;
  }
};



vector_t instanceSNP(SimParameters & s)
{
  
  // Return:
  // 0 Population allele frequency (disease variant)
  // 1 Population allele frequency (marker)

  // 2 Case AA (marker)
  // 3 Case AB (marker)
  
  // 4 Control AA (marker)
  // 5 Control AB (marker)
  
  // 6 P( disease A | marker A ) 
  // 7 P( disease A | marker B )
  

  vector_t freqs(38,0);
  
  // Calculate actual population allele frequency for this SNP
  
  double freq = s.lfreq + CRandom::rand() * ( s.ufreq - s.lfreq ) ;
  
  // And a marker frequency

  double mfreq = par::simul_tags ? 
    s.lmarker + CRandom::rand() * ( s.umarker - s.lmarker ) :
    freq;


  // Handle LD

  double dmax;
  double ld = 0;
  double h11, h12, h21, h22;

  dmax = freq * (1-mfreq);
  if ( (1-freq) * mfreq < dmax) dmax = (1-freq) * mfreq;
  ld = s.dprime * dmax;


  // Haplotype frequencies in general population

  h11 = freq * mfreq + ld;
  h12 = freq * ( 1 - mfreq ) - ld;
  h21 = ( 1 - freq ) * mfreq - ld;
  h22 = ( 1 - freq ) * ( 1 - mfreq ) + ld;


  // Get case and control allele frequencies given GRR, 
  // population frequency and disease frequency
  
  // Need to model 

  //     AM / AM     
  // 2 * AM / Am
  //     Am / Am

  // 2*  AM / aM
  // 4*  AM / am
  // 2*  Am / am

  //     aM / aM
  // 2 * aM / am
  //     am / am
  

  // P( disease | causal variant )

  double f0, f1, f2;
  double g0 = freq * freq;
  double g1 = 2 * freq * ( 1-freq );
  double g2 = 1 - g0 - g1;

  f2 = par::simul_prevalence / ( g0 * s.homOdds + g1 * s.hetOdds + g2 );
  f0 = f2 * s.homOdds;
  f1 = f2 * s.hetOdds;
  

  // P ( disease | diplotype )
  
  double gh_11_11 = h11*h11;
  double gh_11_12 = h11*h12;
  double gh_12_11 = h12*h11;
  double gh_12_12 = h12*h12;

  double gh_11_21 = h11*h21;
  double gh_11_22 = h11*h22;
  double gh_12_21 = h12*h21;
  double gh_12_22 = h12*h22;

  double gh_21_11 = h21*h11;
  double gh_21_12 = h21*h12;
  double gh_22_11 = h22*h11;
  double gh_22_12 = h22*h12;

  double gh_21_21 = h21*h21;
  double gh_21_22 = h21*h22;
  double gh_22_21 = h22*h21;
  double gh_22_22 = h22*h22;

  
  // P( disease | diplotype )

  double fh_11_11 = f0;
  double fh_11_12 = f0;
  double fh_12_11 = f0;
  double fh_12_12 = f0;

  double fh_11_21 = f1;
  double fh_11_22 = f1;
  double fh_12_21 = f1;
  double fh_12_22 = f1;

  double fh_21_11 = f1;
  double fh_21_12 = f1;
  double fh_22_11 = f1;
  double fh_22_12 = f1;

  double fh_21_21 = f2;
  double fh_21_22 = f2;
  double fh_22_21 = f2;
  double fh_22_22 = f2; 
  
  
  double mg0 = mfreq * mfreq;
  double mg1 = 2 * mfreq * ( 1 - mfreq);
  double mg2 = 1 - mg0 - mg1;

  double mf0 = ( f0 * h11*h11
 		 + f1 * 2 * h11 * h21 
		 + f2 * h21*h21 ) / mg0;
  
  double mf1 = ( f0 * 2*h11*h12
		 + f1 * ( 2 * h11 * h22 + 2 * h12 * h21 ) 
		 + f2 * 2*h21*h22 ) / mg1;
  
  double mf2 = ( f0 * h12*h12
		 + f1 * 2 * h12 * h22 
		 + f2 * h22*h22 ) / mg2;
  
    // P(G|X)

  double d0 = mg0 * mf0;
  double d1 = mg1 * mf1;
  double d2 = mg2 * mf2;
  double dSb = d0 + d1 + d2;
  d0 /= dSb;
  d1 /= dSb;
  d2 /= dSb;

  double u0 = mg0 * (1-mf0);
  double u1 = mg1 * (1-mf1);
  double u2 = mg2 * (1-mf2);
  double uSb = u0 + u1 + u2;
  u0 /= uSb;
  u1 /= uSb;
  u2 /= uSb;

  
  // P( diplotype | affected )
  
  double ah_11_11 = fh_11_11 * gh_11_11;
  double ah_11_12 = fh_11_12 * gh_12_11;
  double ah_12_11 = fh_12_11 * gh_12_11;
  double ah_12_12 = fh_12_12 * gh_12_12;
	 	    	       		 	 	 
  double ah_11_21 = fh_11_21 * gh_11_21;
  double ah_11_22 = fh_11_22 * gh_11_22;
  double ah_12_21 = fh_12_21 * gh_12_21;
  double ah_12_22 = fh_12_22 * gh_12_22;
	 	    	       		 	 	 
  double ah_21_11 = fh_21_11 * gh_21_11;
  double ah_21_12 = fh_21_12 * gh_21_12;
  double ah_22_11 = fh_22_11 * gh_22_11;
  double ah_22_12 = fh_22_12 * gh_22_12;
	 	    	       		 	 	 
  double ah_21_21 = fh_21_21 * gh_21_21;
  double ah_21_22 = fh_21_22 * gh_21_22;
  double ah_22_21 = fh_22_21 * gh_22_21;
  double ah_22_22 = fh_22_22 * gh_22_22;
    
  double aS = ah_11_11+ah_11_12+ah_12_11+ah_12_12+
    ah_11_21+ah_11_22+ah_12_21+ah_12_22+
    ah_21_11+ah_21_12+ah_22_11+ah_22_12+
    ah_21_21+ah_21_22+ah_22_21+ah_22_22;

  ah_11_11  /= aS;
  ah_11_12  /= aS;
  ah_12_11  /= aS;
  ah_12_12  /= aS;
  	   
  ah_11_21  /= aS;
  ah_11_22  /= aS;
  ah_12_21  /= aS;
  ah_12_22  /= aS;
  	   
  ah_21_11  /= aS;
  ah_21_12  /= aS;
  ah_22_11  /= aS;
  ah_22_12  /= aS;
  	   
  ah_21_21  /= aS;
  ah_21_22  /= aS;
  ah_22_21  /= aS;
  ah_22_22  /= aS;


  // P( diplotype | unaffected )
  
  double uh_11_11 = (1-fh_11_11) * gh_11_11;
  double uh_11_12 = (1-fh_11_12) * gh_12_11;
  double uh_12_11 = (1-fh_12_11) * gh_12_11;
  double uh_12_12 = (1-fh_12_12) * gh_12_12;
	 	    	       		 	 	 
  double uh_11_21 = (1-fh_11_21) * gh_11_21;
  double uh_11_22 = (1-fh_11_22) * gh_11_22;
  double uh_12_21 = (1-fh_12_21) * gh_12_21;
  double uh_12_22 = (1-fh_12_22) * gh_12_22;
	 	    	       		 	 	 
  double uh_21_11 = (1-fh_21_11) * gh_21_11;
  double uh_21_12 = (1-fh_21_12) * gh_21_12;
  double uh_22_11 = (1-fh_22_11) * gh_22_11;
  double uh_22_12 = (1-fh_22_12) * gh_22_12;
	 	    	       		 	 	 
  double uh_21_21 = (1-fh_21_21) * gh_21_21;
  double uh_21_22 = (1-fh_21_22) * gh_21_22;
  double uh_22_21 = (1-fh_22_21) * gh_22_21;
  double uh_22_22 = (1-fh_22_22) * gh_22_22;
    
  double uS = uh_11_11+uh_11_12+uh_12_11+uh_12_12+
    uh_11_21+uh_11_22+uh_12_21+uh_12_22+
    uh_21_11+uh_21_12+uh_22_11+uh_22_12+
    uh_21_21+uh_21_22+uh_22_21+uh_22_22;

  uh_11_11  /= uS;
  uh_11_12  /= uS;
  uh_12_11  /= uS;
  uh_12_12  /= uS;
  
  uh_11_21  /= uS;
  uh_11_22  /= uS;
  uh_12_21  /= uS;
  uh_12_22  /= uS;
  
  uh_21_11  /= uS;
  uh_21_12  /= uS;
  uh_22_11  /= uS;
  uh_22_12  /= uS;
  
  uh_21_21  /= uS;
  uh_21_22  /= uS;
  uh_22_21  /= uS;
  uh_22_22 /= uS;



  // Return vector

  freqs[0] = freq;   // P(variant)
  freqs[1] = mfreq;  // P(marker)

  freqs[2] = d0;     // P(marker-het|affected)
  freqs[3] = d1;     // P(marker-hom|affected)

  freqs[4] = u0;     // P(marker-het|unaffected)
  freqs[5] = u1;     // P(marker-hom|unaffected)

  
  freqs[6] =     ah_11_11 ;
  freqs[7] =     ah_11_12 ;
  freqs[8] =     ah_12_11 ;
  freqs[9] =     ah_12_12 ;
  
  freqs[10] =    ah_11_21 ;
  freqs[11] =    ah_11_22 ;
  freqs[12] =    ah_12_21 ;
  freqs[13] =    ah_12_22 ;
  
  freqs[14] =    ah_21_11 ;
  freqs[15] =    ah_21_12 ;
  freqs[16] =    ah_22_11 ;
  freqs[17] =    ah_22_12 ;
  
  freqs[18] =    ah_21_21 ;
  freqs[19] =    ah_21_22 ;
  freqs[20] =    ah_22_21 ;
  freqs[21] =    ah_22_22 ;

  // Control diplotype freqs
  freqs[22] =     uh_11_11 ;
  freqs[23] =     uh_11_12 ;
  freqs[24] =     uh_12_11 ;
  freqs[25] =     uh_12_12 ;
  
  freqs[26] =    uh_11_21 ;
  freqs[27] =    uh_11_22 ;
  freqs[28] =    uh_12_21 ;
  freqs[29] =    uh_12_22 ;
  
  freqs[30] =    uh_21_11 ;
  freqs[31] =    uh_21_12 ;
  freqs[32] =    uh_22_11 ;
  freqs[33] =    uh_22_12 ;
  
  freqs[34] =    uh_21_21 ;
  freqs[35] =    uh_21_22 ;
  freqs[36] =    uh_22_21 ;
  freqs[37] =    uh_22_22 ;

  return freqs;
}

vector_t instanceSNP_QT(SimParametersQT & s)
{
  
  // Return:

  // 0 Population allele frequency (disease variant)
  // 1 Population allele frequency (marker)
  
  // 16 haplotype frequencies
  

  vector_t freqs(18,0);
  
  // Calculate actual population allele frequency for this SNP
  
  double freq = s.lfreq + CRandom::rand() * ( s.ufreq - s.lfreq ) ;
  
  // And a marker frequency

  double mfreq = par::simul_tags ? 
    s.lmarker + CRandom::rand() * ( s.umarker - s.lmarker ) :
    freq;

  // Given the specified allele frequency, now calcuate the 
  // additive genetic value and domiance deviation

  double p = freq;
  double q = 1-p;
  double a = sqrt( ( s.variance ) 
		   / ( (2*p*q)* (1+s.dom*(q-p))*(1+s.dom*(q-p)) + (2*p*q*s.dom)*(2*p*q*s.dom) ) );
  double d = s.dom * a;
 
  // Mean center

 
  s.gBB =  a  -(a*(p-(1-p))+ (2*p*(1-p)*d));
  s.gAB =  d  -(a*(p-(1-p))+ (2*p*(1-p)*d));
  s.gAA = -a  -(a*(p-(1-p))+ (2*p*(1-p)*d));
  
//   cout << "p = " << p << "\n";
//   cout <<"G = " << s.gBB << " " 
//        << s.gAB << " " 
//       << s.gAA << "\n";

    


  // Handle LD

  double dmax;
  double ld = 0;
  double h11, h12, h21, h22;

  dmax = freq * (1-mfreq);
  if ( (1-freq) * mfreq < dmax) dmax = (1-freq) * mfreq;
  ld = s.dprime * dmax;


  // Haplotype frequencies in general population

  h11 = freq * mfreq + ld;
  h12 = freq * ( 1 - mfreq ) - ld;
  h21 = ( 1 - freq ) * mfreq - ld;
  h22 = ( 1 - freq ) * ( 1 - mfreq ) + ld;

  double h_11_11 = h11*h11;
  double h_11_12 = h11*h12;
  double h_12_11 = h12*h11;
  double h_12_12 = h12*h12;

  double h_11_21 = h11*h21;
  double h_11_22 = h11*h22;
  double h_12_21 = h12*h21;
  double h_12_22 = h12*h22;

  double h_21_11 = h21*h11;
  double h_21_12 = h21*h12;
  double h_22_11 = h22*h11;
  double h_22_12 = h22*h12;

  double h_21_21 = h21*h21;
  double h_21_22 = h21*h22;
  double h_22_21 = h22*h21;
  double h_22_22 = h22*h22;


  // Return vector

  freqs[0] = freq;   // P(variant)
  freqs[1] = mfreq;  // P(marker)
  
  freqs[2] =     h_11_11 ;
  freqs[3] =     h_11_12 ;
  freqs[4] =     h_12_11 ;
  freqs[5] =     h_12_12 ;
  
  freqs[6] =    h_11_21 ;
  freqs[7] =    h_11_22 ;
  freqs[8] =    h_12_21 ;
  freqs[9] =    h_12_22 ;
  
  freqs[10] =    h_21_11 ;
  freqs[11] =    h_21_12 ;
  freqs[12] =    h_22_11 ;
  freqs[13] =    h_22_12 ;
  
  freqs[14] =    h_21_21 ;
  freqs[15] =    h_21_22 ;
  freqs[16] =    h_22_21 ;
  freqs[17] =    h_22_22 ;

  return freqs;
}

void Plink::simulateSNPs()
{

  // Read in SNP parameters

  // Number of SNPs
  // Lower allele frequency for '1' (versus '2') allele
  // Upper allele frequency (population)
  // Odds ratio ('1' allele)

  checkFileExists(par::simul_file);
  printLOG("Reading simulation parameters from [ " 
	   + par::simul_file + " ]\n");
  
  printLOG("Writing SNP population frequencies to [ " 
	   + par::output_file_name + ".simfreq ]\n");
  
  ofstream SOUT( ( par::output_file_name+".simfreq").c_str(), ios::out);
  
  ifstream SIM;
  SIM.open( par::simul_file.c_str(), ios::in );
  
  if ( par::simul_label != "" ) 
    par::simul_label += "-";

  vector<SimParameters> sp;

  while ( ! SIM.eof() )
    {

      SimParameters s;
      
      vector<string> tokens = tokenizeLine( SIM );
      
      if ( tokens.size() == 0 )
	continue;
      
      if ( par::simul_tags )
	{
	  	  
	  if ( tokens.size() != 9 )
	    error("Problem with format of simulation parameter file: expecting 9 fields\n");
	  
	  if( ! from_string<int>(s.nsnp , tokens[0] , std::dec) )
	    error("Expecting numeric value for 1st field, # SNPs\n");
	  
	  s.name = tokens[1];

	  if( ! from_string<double>(s.lfreq , tokens[2] , std::dec) )
	    error("Expecting numeric value for 3rd field, lower variant freq.\n");

	  if( ! from_string<double>(s.ufreq , tokens[3] , std::dec) )
	    error("Expecting numeric value for 4th field, upper variant freq.\n");
	  
	  if( ! from_string<double>(s.lmarker , tokens[4] , std::dec) )
	    error("Expecting numeric value for 5th field, lower marker freq.\n");

	  if( ! from_string<double>(s.umarker , tokens[5] , std::dec) )
	    error("Expecting numeric value for 6th field, upper marker freq.\n");

	  if( ! from_string<double>(s.dprime , tokens[6] , std::dec) )
	    error("Expecting numeric value for 7th field, d-prime\n");

	  if( ! from_string<double>(s.hetOdds , tokens[7] , std::dec) )
	    error("Expecting numeric value for 8th field, het odds\n");
	  
	  if (  ! from_string<double>( s.homOdds , tokens[8] , std::dec ) )
	    s.homOdds = s.hetOdds * s.hetOdds;
	  
	}
      else
	{      
	  if ( tokens.size() != 6 )
	    error("Problem with format of simulation parameter file: expecting 6 fields\n");

	  if( ! from_string<int>(s.nsnp , tokens[0] , std::dec) )
	    error("Expecting numeric value for first field, # SNPs\n");
	  
	  s.name = tokens[1];
	  
	  if( ! from_string<double>(s.lfreq , tokens[2] , std::dec) )
	    error("Expecting numeric value for 3rd field, lower variant freq.\n");
	  
	  if( ! from_string<double>(s.ufreq , tokens[3] , std::dec) )
	    error("Expecting numeric value for 4th field, upper variant freq.\n");

	  s.lmarker = s.lfreq;
	  s.umarker = s.ufreq;
	  s.dprime = 1;
	  
	  if( ! from_string<double>(s.hetOdds , tokens[4] , std::dec) )
	    error("Expecting numeric value for 5th field, het odds\n");
	  
	  if (  ! from_string<double>( s.homOdds , tokens[5] , std::dec ) )
	    s.homOdds = s.hetOdds * s.hetOdds;
	  
	}

      // Read odds ratio; unless a specific number is given, assume
      // multiplicative

            
      sp.push_back(s);           
      
      if ( SIM.eof() )
	break;
      
    }
  
  SIM.close();



  ////////////////////////////////////////////
  // Make room for total number of SNPs, etc
  
  int tsnp = 0;
  
  for (int s=0; s<sp.size(); s++)
    tsnp += sp[s].nsnp;
  
  if ( par::simul_haps )
    tsnp *= 2;

  int nind = par::simul_ncases + par::simul_ncontrols;
      
  printLOG("Read " + int2str( sp.size() ) 
	   + " sets of SNPs, specifying " 
	   + int2str(tsnp) 
	   + " SNPs in total\n");
  printLOG("Simulating " 
	   + int2str( par::simul_ncases ) 
	   + " cases and " 
	   + int2str( par::simul_ncontrols ) 
	   + " controls\n");
  printLOG("Assuming a disease prevalence of " 
	   + dbl2str( par::simul_prevalence ) + "\n");
  

  if ( par::simul_haps )
    printLOG("Simulating causal variants and markers\n\n");
  else if ( par::simul_tags )
    printLOG("Simulating markers in indirect association\n\n");
  else
    printLOG("Simulating disease variants (direct association)\n\n");
  
  vector<vector<bool> > haps;
  for (int cv1=0; cv1<2; cv1++)
    for (int cv2=0; cv2<2; cv2++)
      for (int mk1=0; mk1<2; mk1++)
	for (int mk2=0; mk2<2; mk2++)
	  {
	    vector<bool> t(4);
	    t[0] = cv1;
	    t[1] = mk1;
	    t[2] = cv2;
	    t[3] = mk2;
	    haps.push_back(t);	    
	  }
  
  int pos = 0;
  
  for (int s=0; s<sp.size(); s++)  
    for (int l=0; l < sp[s].nsnp; l++ )
      {
	
	// Causal variant

	Locus * loc = new Locus;
	
	// Optionally add
	if ( sp[s].nsnp > 1 ) 
	  loc->name = sp[s].name+"_"+int2str(l);
	else
	  loc->name = sp[s].name;
		
	loc->chr = 1;
	loc->allele1 = "D";
	loc->allele2 = "d";
	loc->bp = ++pos; 
	loc->pos = 0;
	
	locus.push_back(loc);
	
	CSNP * newset = new CSNP;
	
	newset->one.resize(nind);
	newset->two.resize(nind);

	if ( par::simul_haps )
	  {
	    Locus * loc2 = new Locus;	
	    loc2->name = loc->name + "_M";
	    loc2->chr = 1;
	    loc2->allele1 = "A";
	    loc2->allele2 = "B";
	    loc2->bp = ++pos; 
	    loc2->pos = 0;
	    locus.push_back(loc2);
	  }

	
	  
	// Sample case and control population genotype frequencies
	
	vector_t f = instanceSNP(sp[s]);      
	
	// f Information
	
	// 0 Population allele frequency (disease variant/marker)
	// 1 Population allele frequency (marker)
	// 2 Case AA (marker)
	// 3 Case AB (marker)
	// 4 Control AA (marker)
	// 5 Control AB (marker)
	
	// 6+ full diplotype frequencies, for 
	//    cases and controls

	
	if ( par::simul_tags )
	  {
	    SOUT << 1 << " "
		 << loc->name << "\t"
		 << f[0] << " " << f[0] << "\t" 
		 << f[1] << " " << f[1] << "\t" 
		 << sp[s].dprime << "\t"
		 << sp[s].hetOdds << "\t"
		 << sp[s].homOdds << "\n";
	  }
	else
	  {	  
	    SOUT << 1 << " "
		 << loc->name << "\t"
		 << f[0] << " " << f[0] << "\t" 
		 << sp[s].hetOdds << "\t" 
		 << sp[s].homOdds << "\n";	 
	  }
	
	

	if ( ! par::simul_haps )
	  {
	    
	    // Simulate only a single SNP (either the marker, 
	    // if --simulate-tags, otherwise the CV itself
	    
	    // Genotype frequencies in cases and controls
	
	    const double caseAA = f[2];
	    const double caseAB = f[2] + f[3];
	    
	    const double contAA = f[4];
	    const double contAB = f[4] + f[5];
	    
	    
	    //////////////////////////////////////////////////
	    // Generate each individual, simulating genotypes
	    // rather than alleles
	    
	    for ( int i = 0 ; i < nind ; i++ ) 
	      {
		
		// Simple missingness
		
		if ( CRandom::rand() < sp[s].missing ) 
		  {
		    newset->one[i] = true;
		    newset->two[i] = false;
		  }
		else
		  {
		    
		    bool isCase = i < par::simul_ncases ? true : false;
		    
		    double r = CRandom::rand();
		    
		    int g = 0;
		    		    
		    if ( isCase ) 
		      {
			if ( r > caseAB )
			  g = 2;		  
			else if ( r > caseAA )
			  g = 1;
		      }
		    else
		      {
			if ( r > contAB )
			  g = 2;
			else if ( r > contAA )
			  g = 1;
			
		      }
		    
		    
		    if ( g == 2 ) 
		      {
			newset->one[i] = false;
			newset->two[i] = false;
		      }
		    else if ( g == 1 ) 
		      {
			newset->one[i] = false;
			newset->two[i] = true;
		      }
		    else
		      {
			newset->one[i] = true;
			newset->two[i] = true;
		      }
		    
		  }
	      }
	    
	    SNP.push_back(newset);
	    
	  }
	else
	  {

	    
	    ///////////////////////////////////
	    // Simulate diplotype pair
	    
	    CSNP * newset2 = new CSNP;	    
	    
	    newset2->one.resize(nind);
	    newset2->two.resize(nind);	    
	    
	    // Cases: 6 to 21
	    // Controls: 22 to 37

	    vector_t freqA;
	    vector_t freqU;
	    double cumA = 0, cumU = 0;

	    for (int j=6; j<=21; j++)
	      {
		cumA += f[j];
		freqA.push_back(cumA);
	      }

	    for (int j=22; j<=37; j++)
	      {
		cumU += f[j];
		freqU.push_back(cumU);
	      }

	    
	    //////////////////////////////
	    // Generate each individual, 
	    
	    for ( int i = 0 ; i < nind ; i++ ) 
	      {
		
		// Simple missingness
		bool miss_marker = false;
		bool miss_causal = false;

		if ( CRandom::rand() < sp[s].missing ) 
		  miss_marker = true;

		if ( CRandom::rand() < sp[s].missing ) 
		  miss_causal = true;

		
		bool isCase = i < par::simul_ncases ? true : false;
		
		// Simulate diplotype

		double r = CRandom::rand();
		int h = 0;
		
		for ( int j=14;j>=0;j--)
		  {
		    if ( isCase )
		      {
			if ( r > freqA[j] )
			  {
			    h = j+1;
			    break;
			  }
		      }
		    else 
		      {
			if ( r > freqU[j] )
			  {
			    h = j+1;
			    break;
			  }
		      }
		  }
		
		
		// We now have selected 'h', a number between 0 and 15
		
		vector<bool> & hp = haps[h];
		
		//cout << "h = " << h << "\n";
		// 		    cout << "size haps=" << haps.size() << "\n";
		// 		    display(freqA);
		
		    //////////////////////////
		    // Set both genotypes
		
		    if ( miss_marker )
		      {
			newset->one[i] = true;
			newset->two[i] = false;
		      }
		    else
		      {
			
			if ( hp[0] && hp[2] ) 
			  {
			    newset->one[i] = true;
			    newset->two[i] = true;
			  }
			else if ( (!hp[0]) && (!hp[2]) ) 
			  {
			    newset->one[i] = false;
			    newset->two[i] = false;
			  }
			else
			  {
			    newset->one[i] = false;
			    newset->two[i] = true;
			  }
		      }
		    
		    if ( miss_causal )
		      {
			newset2->one[i] = true;
			newset2->two[i] = false;
		      }
		    else
		      {
			if ( hp[1] && hp[3] ) 
			  {
			    newset2->one[i] = true;
			    newset2->two[i] = true;
			  }
			else if ( (!hp[1]) && (!hp[3]) ) 
			  {
			    newset2->one[i] = false;
			    newset2->two[i] = false;
			  }
			else
			  {
			    newset2->one[i] = false;
			    newset2->two[i] = true;
			  }


		      }
		    

	      }
	      
	    
	    // Add markers and then CV
	    SNP.push_back(newset);
	    SNP.push_back(newset2);

	  }
	
	// Next SNP/SNP-pair to simulate 
      }
  
  
  // Phenotypes
  
  for (int i=0;i<nind;i++)
    {
      Individual * person = new Individual;
      person->fid = person->iid = par::simul_label + "per"+int2str(i);
      person->missing = false;
      person->pat = "0";
      person->mat = "0";

      if ( i < par::simul_ncases ) 
	person->phenotype = 2;
      else
	person->phenotype = 1;

      person->sex = false;
      person->sexcode = "2";
      sample.push_back(person);      
    }
  
  SOUT.close();


}






void Plink::simulateSNPs_QT()
{

  par::qt = true;
  par::bt = false;

  // Read in SNP parameters
  // for QUANTITATIVE TRAIT simulations

  // Number of SNPs  
  // Lower allele frequency for '1' (versus '2') allele
  // Upper allele frequency (population)
  // Additive genetic effect
  // Dominance deviation

  checkFileExists(par::simul_file);

  printLOG("Reading QT simulation parameters from [ " 
	   + par::simul_file + " ]\n");
  
  printLOG("Writing SNP population frequencies to [ " 
	   + par::output_file_name + ".simfreq ]\n");
  
  ofstream SOUT( ( par::output_file_name+".simfreq").c_str(), ios::out);
  
  ifstream SIM;
  SIM.open( par::simul_file.c_str(), ios::in );
  
  if ( par::simul_label != "" ) 
    par::simul_label += "-";

  vector<SimParametersQT> sp;

  double totvar = 0;

  while ( ! SIM.eof() )
    {

      SimParametersQT s;
      
      vector<string> tokens = tokenizeLine( SIM );
      
      if ( tokens.size() == 0 )
	continue;
      
      
      if ( par::simul_tags )
	{
	  
	  if ( tokens.size() != 9 )
	    error("Problem with format of simulation parameter file: expecting 9 fields\n");
	  
	  if( ! from_string<int>(s.nsnp , tokens[0] , std::dec) )
	    error("Expecting numeric value for 1st field, # SNPs\n");
	  
	  s.name = tokens[1];

	  if( ! from_string<double>(s.lfreq , tokens[2] , std::dec) )
	    error("Expecting numeric value for 3rd field, lower variant freq.\n");

	  if( ! from_string<double>(s.ufreq , tokens[3] , std::dec) )
	    error("Expecting numeric value for 4th field, upper variant freq.\n");
	  
	  if( ! from_string<double>(s.lmarker , tokens[4] , std::dec) )
	    error("Expecting numeric value for 5th field, lower marker freq.\n");

	  if( ! from_string<double>(s.umarker , tokens[5] , std::dec) )
	    error("Expecting numeric value for 6th field, upper marker freq.\n");

	  if( ! from_string<double>(s.dprime , tokens[6] , std::dec) )
	    error("Expecting numeric value for 7th field, d-prime\n");

	  if( ! from_string<double>(s.variance , tokens[7] , std::dec) )
	    error("Expecting numeric value for 8th field, variance\n");
	  
	  if (  ! from_string<double>( s.dom , tokens[8] , std::dec ) )
	    error("Expecting numeric value for 9th field, dom \n");
	  
	}
      else
	{      


	  if ( tokens.size() != 6 )
	    error("Problem with format of simulation parameter file: expecting 6 fields\n");
      
	  if( ! from_string<int>(s.nsnp , tokens[0] , std::dec) )
	    error("Expecting numeric value for first field, # SNPs\n");
	  
	  s.name = tokens[1];
	  
	  if( ! from_string<double>(s.lfreq , tokens[2] , std::dec) )
	    error("Expecting numeric value for 3rd field, lower variant freq.\n");
	  
	  if( ! from_string<double>(s.ufreq , tokens[3] , std::dec) )
	    error("Expecting numeric value for 4th field, upper variant freq.\n");

	  if( ! from_string<double>(s.variance , tokens[4] , std::dec) )
	    error("Expecting numeric value for 5th field, variance\n");
      
	  if (  ! from_string<double>( s.dom , tokens[5] , std::dec ) )
	    error("Expecting numeric value for 6th field, dom \n");
	  
	  s.lmarker = s.lfreq;
	  s.umarker = s.ufreq;
	  s.dprime = 1;
	  
	}
      
      // Keep track of total QTL variance
      totvar += s.nsnp * s.variance;

      sp.push_back(s);           
      
      if ( SIM.eof() )
	break;
      
    }
  
  SIM.close();
  
  
  
  ////////////////////////////////////////////
  // Make room for total number of SNPs, etc
  
  int tsnp = 0;
  
  for (int s=0; s<sp.size(); s++)
    tsnp += sp[s].nsnp;
  
  if ( par::simul_haps )
    tsnp *= 2;
  
  int nind = par::simul_ncases;
      
  printLOG("Read " + int2str( sp.size() ) 
	   + " sets of SNPs, specifying " 
	   + int2str(tsnp) 
	   + " SNPs in total\n");
  printLOG("Simulating QT for " 
	   + int2str( par::simul_ncases ) 
	   + " individuals \n");
  printLOG("Total QTL variance is " 
	   + dbl2str( totvar ) + "\n");

  if ( par::simul_haps )
    printLOG("Simulating causal variants and markers\n\n");
  else if ( par::simul_tags )
    printLOG("Simulating markers in indirect association\n\n");
  else
    printLOG("Simulating disease variants (direct association)\n\n");
  
  if ( totvar > 1 )
    error("Specific QTL variance is greater than 100%");
  
  // Phenotypes
  
  for (int i=0;i<nind;i++)
    {
      Individual * person = new Individual;
      person->fid = person->iid = par::simul_label + "per"+int2str(i);
      person->missing = false;
      person->pat = "0";
      person->mat = "0";
      person->sex = false;
      person->sexcode = "2";

      // Residual variance component
      person->phenotype = rnorm() * sqrt( 1 - totvar );

      sample.push_back(person);      
    }

  
  // Genotypes
  
  vector<vector<bool> > haps;
  for (int cv1=0; cv1<2; cv1++)
    for (int cv2=0; cv2<2; cv2++)
      for (int mk1=0; mk1<2; mk1++)
	for (int mk2=0; mk2<2; mk2++)
	  {
	    vector<bool> t(4);
	    t[0] = cv1;
	    t[1] = mk1;
	    t[2] = cv2;
	    t[3] = mk2;
	    haps.push_back(t);	    
	  }
  
  int pos = 0;
  
  for (int s=0; s<sp.size(); s++)  
    for (int l=0; l < sp[s].nsnp; l++ )
      {
	
	// Causal variant

	Locus * loc = new Locus;
	
	// Optionally add
	if ( sp[s].nsnp > 1 ) 
	  loc->name = sp[s].name+"_"+int2str(l);
	else
	  loc->name = sp[s].name;
		
	loc->chr = 1;
	loc->allele1 = "H";
	loc->allele2 = "L";
	loc->bp = ++pos; 
	loc->pos = 0;
	
	locus.push_back(loc);
	
	CSNP * newset = new CSNP;
	
	newset->one.resize(nind);
	newset->two.resize(nind);

	if ( par::simul_haps )
	  {
	    Locus * loc2 = new Locus;	
	    loc2->name = loc->name + "_M";
	    loc2->chr = 1;
	    loc2->allele1 = "A";
	    loc2->allele2 = "B";
	    loc2->bp = ++pos; 
	    loc2->pos = 0;
	    locus.push_back(loc2);
	  }

	
	  
	// Get haplotype frequencies
	
	vector_t f = instanceSNP_QT(sp[s]);      
	
	// f Information
	
	// 0 Population allele frequency (disease variant/marker)
	// 1 Population allele frequency (marker)
	// 2+  16 haplotype frequencies
	
	if ( par::simul_tags )
	  {
	    SOUT << 1 << " "
		 << loc->name << "\t"
		 << f[0] << " " << f[0] << "\t" 
		 << f[1] << " " << f[1] << "\t" 
		 << sp[s].dprime << "\t"
		 << sp[s].variance << "\t"
		 << sp[s].dom << "\n";
	  }
	else
	  {	  
	    SOUT << 1 << " "
		 << loc->name << "\t"
		 << f[0] << " " << f[0] << "\t" 
		 << sp[s].variance << "\t" 
		 << sp[s].dom << "\n";	 
	  }
	
	

	if ( ! par::simul_haps )
	  {
	    
	    // Simulate only a single SNP (either the marker, 
	    // if --simulate-tags, otherwise the CV itself
	    
	    // Genotype frequencies in cases and controls
	
	    const double freq  = f[0];
	    
	    
	    //////////////////////////////////////////////////
	    // Generate each individual, simulating genotypes
	    // rather than alleles
	    
	    for ( int i = 0 ; i < nind ; i++ ) 
	      {
		
		// Simple missingness
		
		if ( CRandom::rand() < sp[s].missing ) 
		  {
		    newset->one[i] = true;
		    newset->two[i] = false;
		  }
		else
		  {
		    		    		    
		    int g = 0;
		    if ( CRandom::rand() > freq ) ++g;
		    if ( CRandom::rand() > freq ) ++g;
		    
		    if ( g == 2 ) 
		      {
			newset->one[i] = false;
			newset->two[i] = false;
			sample[i]->phenotype += sp[s].gAA;
		      }
		    else if ( g == 1 ) 
		      {
			newset->one[i] = false;
			newset->two[i] = true;
			sample[i]->phenotype += sp[s].gAB;
		      }
		    else
		      {
			newset->one[i] = true;
			newset->two[i] = true;
			sample[i]->phenotype += sp[s].gBB;
			
		      }
		    
		  }
	      }
	    
	    SNP.push_back(newset);
	    
	  }
	else
	  {

	    
	    ///////////////////////////////////
	    // Simulate diplotype pair
	    
	    CSNP * newset2 = new CSNP;	    
	    
	    newset2->one.resize(nind);
	    newset2->two.resize(nind);	    
	    
	    // Cases: 2 to 17

	    
	    vector_t freqA;
	    vector_t freqU;
	    double cumA = 0, cumU = 0;

	    for (int j=2; j<=17; j++)
	      {
		cumA += f[j];
		freqA.push_back(cumA);
	      }

	    
	    //////////////////////////////
	    // Generate each individual, 
	    
	    for ( int i = 0 ; i < nind ; i++ ) 
	      {
		
		// Simple missingness
		bool miss_marker = false;
		bool miss_causal = false;

		if ( CRandom::rand() < sp[s].missing ) 
		  miss_marker = true;

		if ( CRandom::rand() < sp[s].missing ) 
		  miss_causal = true;

		
		// Simulate diplotype

		double r = CRandom::rand();
		int h = 0;
		
		for ( int j=14;j>=0;j--)
		  {
		    if ( r > freqA[j] )
		      {
			h = j+1;
			break;
		      }
		  }
		
		
		// We now have selected 'h', a number between 0 and 15
		
		vector<bool> & hp = haps[h];
		
		cout << "h = " << h << "\n";
		cout << "size haps=" << haps.size() << "\n";
		display(freqA);
		
		//////////////////////////
		// Set both genotypes
		
		if ( miss_marker )
		  {
		    newset->one[i] = true;
		    newset->two[i] = false;
		  }
		else
		  {
		    
		    if ( hp[0] && hp[2] ) 
		      {
			newset->one[i] = true;
			newset->two[i] = true;
		      }
		    else if ( (!hp[0]) && (!hp[2]) ) 
		      {
			newset->one[i] = false;
			newset->two[i] = false;
		      }
		    else
		      {
			newset->one[i] = false;
			newset->two[i] = true;
		      }
		  }
		
		if ( miss_causal )
		  {
		    newset2->one[i] = true;
		    newset2->two[i] = false;
		  }
		else
		  {
		    if ( hp[1] && hp[3] ) 
		      {
			newset2->one[i] = true;
			newset2->two[i] = true;
			sample[i]->phenotype += sp[s].gAA;
		      }
		    else if ( (!hp[1]) && (!hp[3]) ) 
		      {
			newset2->one[i] = false;
			newset2->two[i] = false;
			sample[i]->phenotype += sp[s].gAB;
		      }
		    else
		      {
			newset2->one[i] = false;
			newset2->two[i] = true;
			sample[i]->phenotype += sp[s].gBB;
		      }
		    
		  }
		
	      }
	      
	    
	    // Add markers and then CV
	    SNP.push_back(newset);
	    SNP.push_back(newset2);

	  }
	
	// Next SNP/SNP-pair to simulate 
      }
  
    
  SOUT.close();


}







#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/fft.h"

#include "glsol.hpp"
#include "matep.hpp"


void glsol::case_8() {
   
    pi = 0; 
    deltaPi =0;
    djAaj = 0;   
    phaseMarker = 0.;
    
    // set all sites to be normal phase with thermal noise
    A[ALL] = sqrt(config.variance_sigma) * A[X].gaussian_random();

    real_t Tcp_mK = MP.Tcp_mK(config.Inip);
    
    onsites(ALL) {

      matep::Matep MPonsites;
      if (T[X] < Tcp_mK)
	{
         foralldir(al) foralldir(i)
	   {	
	    if ((al==0) && (i==0))
	      { A[X].e(al,i).re = 1.; }
	    else if ((al==0) && (i==1))
	      { A[X].e(al,i).im = 1.; } // put bulk A-phase elements into OP
           } // doralldir end here

         A[X]=A[X] * (MPonsites.gap_A_td(config.Inip, T[X])/sqrt(2.));
	 // hila::out0 << "A[x] is " << A[X].e(0,0).re << std::endl;
	  
       } // Temeprature judgement block
      
    } // onsites(ALL) end here
    hila::out0 << "kBTC/f0p is " << MP.kBTCf0p_ratio(config.Inip)
               << "\n" 
	       << "Thermal fluctuation energy kB*TC is " << MP.kBTC(config.Inip) << "J."
               << "\n"      
	       << "1/3 N(0) (xi0GL)^3 (kB TC)^2 is " << MP.f0p(config.Inip) << "J."
               << "\n"      
	       << "xi0GL is " << MP.xi0GLp(config.Inip) << "m."
               << "\n"      
	       << "N(0) is " << MP.N0p(config.Inip) << "J^-1. m^-3."
               << "\n"      
	       << "mEff is " << MP.mEffp(config.Inip) << "kg."
               << "\n"      
	       << "Fermi velocity vFp is "<< MP.vFp(config.Inip) << "m.s^-1"
               << "\n"      
             // << ", MP.hbar() is " << MP.hbar()
             // << ", MP.hbar3() is " << MP.hbar3()
	       << std::endl;

    hila::out0 << "OP field initialized according to the hotblob profile! " << std::endl;

} // case_8() call end here


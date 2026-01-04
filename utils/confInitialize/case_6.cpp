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


void glsol::case_6() {
       
    pi = 0;
    deltaPi =0;
    djAaj = 0;
    phaseMarker = 0;
        
    real_t gapA = MP.gap_A_td(config.Inip, config.IniT);
    real_t gapB = MP.gap_B_td(config.Inip, config.IniT);
    real_t tb = config.Inilc;//config.IniT/ MP.Tcp_mK(config.Inip);
    hila::out0 << "Gap A: " << gapA <<"\n";
    hila::out0 << "Gap B: " << gapB <<"\n";
    onsites(ALL) {
      real_t d=sqrt(pow(X.coordinates()[0]-config.lx/2.0,2.0)+pow(X.coordinates()[1]-config.ly/2.0,2.0)+pow(X.coordinates()[2]-config.lz/2.0,2.0));
      if(d<tb){
        foralldir(al) foralldir(i){
          A[X].e(al,i) = sqrt(config.IniMod) * hila::gaussian_random<Complex<real_t>>();
        }
      }
      else{
	//A[X].gaussian_random(config.IniMod);
      
	foralldir(al) foralldir(i){
	  if ((al==0) && (i==0)) {
	    A[X].e(al,i).re = 1.;
	  }
	  else if ((al==0) && (i==1)) {
	    A[X].e(al,i).im = 1.;
	  }
	}
	A[X] = gapA * A[X]/sqrt(2.0);
      }
      
    }
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
    
    hila::out0 << "Aphase_partial is created \n";

} // case_6() call end here


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


void glsol::case_5() {
   
    pi = 0;
    deltaPi =0;
    djAaj = 0;
    phaseMarker = 0.F;
    
    real_t gap = MP.gap_B_td(config.Inip, config.IniT);
    hila::out0 << " Gap B : " << gap << std::endl;
    onsites(ALL) {
      foralldir(al)foralldir(i){
	A[X].e(al,i) = sqrt(config.variance_sigma) * hila::gaussian_random<Complex<real_t>>();	

	if (al==i){
	  A[X].e(al,i).re = 1.0 + A[X].e(al,i).re; 
	  A[X].e(al,i).im = 0.0 + A[X].e(al,i).im;
	}	
      }
      A[X] = (gap/sqrt(3.)) * A[X];
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
    
    hila::out0 << "Pure B phase \n";

} // case_5() call end here


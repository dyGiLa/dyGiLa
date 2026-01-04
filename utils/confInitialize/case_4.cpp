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

void glsol::case_4() {
   
    pi = 0;
    deltaPi =0;
    djAaj = 0;    
    phaseMarker = 0.F;
    
    hila::out0 << "gapA = " << MP.gap_A_td(config.Inip, config.IniT) << "at p = " << config.Inip << ", T = " << config.IniT
               << "\n"
               << "gapB = " << MP.gap_B_td(config.Inip, config.IniT) << "at p = " << config.Inip << ", T = " << config.IniT
               << std::endl;
    onsites(ALL) {
      //hila::out0 << "this is case 4" << std::endl;
      foralldir(al) foralldir(i){
	A[X].e(al,i) = sqrt(config.variance_sigma) * hila::gaussian_random<Complex<real_t>>();
      } // doralldir end here
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
  
    hila::out0 << " normal-phase-complex created! " << std::endl;;
       
} // case_4() call end here


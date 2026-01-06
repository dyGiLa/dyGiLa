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


void glsol::case_10() {
   
    pi = 0;
    deltaPi =0;
    djAaj = 0;
    phaseMarker = 0.F;
    
    real_t gap = MP.gap_B_td(config.Inip, config.IniT);
    hila::out0 << " Gap B : " << gap << std::endl;
    onsites(ALL)
      {
	// do the coordinate transformation
	auto x = X.coordinate(e_x) - config.lx/2.f;

	// stripe B half period
	auto shp = config.lx/4.f;

	if (
	    ((x >= 0.f) && (x-shp <= 0.f))
	    || ((x < 0.f) && (x-(-shp) >= 0.f))
	   )
	  {
	    foralldir(al)foralldir(i){
	      A[X].e(al,i) = sqrt(config.variance_sigma) * hila::gaussian_random<Complex<real_t>>();
	      if (al==i){
	        A[X].e(al,i).re += 1.0; 
	        A[X].e(al,i).im += 0.0 ;
	      }	
            }
            A[X] = (gap/sqrt(3.)) * A[X];
	  } // {1,1,1} B-phase
	else if (
	         ((x > 0.f) && (x-shp > 0.f))
	         || ((x < 0.f) && (x-(-shp) < 0.f))
                )
	  {
	    foralldir(al)foralldir(i){
	      A[X].e(al,i) = sqrt(config.variance_sigma) * hila::gaussian_random<Complex<real_t>>();
	      if (al==i && al!=2){
	        A[X].e(al,i).re += 1.0; 
	        A[X].e(al,i).im += 0.0 ;
	      }
	      else if ((al==i && al==2)) {
	        A[X].e(al,i).re += -1.0; 
	        A[X].e(al,i).im += 0.0 ;
	      }		
            }
            A[X] = (gap/sqrt(3.)) * A[X];
	  } // {1, 1, -1} B-phase
	
      } // onsites block ends here

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
    
    hila::out0 << "Wiman2016 Stripe B-phase conf is intialized! " << std::endl;

} // case_10() call end here


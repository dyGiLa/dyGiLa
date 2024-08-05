#define _USE_MATH_DEFINES
#define USE_ASCENT 
#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
//#include <math.h>
#include <assert.h>

#include "plumbing/hila.h"
#include "plumbing/fft.h"

#include "glsol.hpp"
#include "matep.hpp"


void glsol::initialize() {

  Matep MP;
   
  int Nx = config.lx;
  int Ny = config.ly;
  int Nz = config.lz;
  
  real_t dx = config.dx;

  real_t ttc = MP.Tcp_mK(config.Inip);
  hila::out0 <<"T_AB: "<<MP.tAB_RWS(config.Inip)*ttc<<"\n";
  
  switch (config.initialCondition) {
    
  case 0: {
    pi = 0;                            
    real_t gap = MP.gap_B_td(config.Inip, config.IniT);
    onsites(ALL) {                     
      A[X] = hila::gaussrand();
      A[X] = gap * A[X]/A[X].norm();   
    }

    hila::out0 << "Components randomly created \n";

    break;
    }
  case 1: {
    auto kA = A;
    real_t gap = config.IniMod;        //MP.gap_B_td(Tp[1], Tp[0]);
    real_t lc = config.Inilc;          //1.0/sqrt(abs(config.alpha));

    hila::out0 << "Correlation length in ICs: "<< lc <<"\n";
    
    onsites (ALL) {
            real_t constant = pow(gap, 2.0) * pow(2.0 * M_PI, 1.5) * pow(lc, 3.0)/(9.0 * Nx * Ny * Nz * dx * dx * dx);
            real_t kSqu;
            real_t std;
            kSqu = 0.0;
            auto k = X.coordinates();

            kSqu = pow(sin(M_PI * k.e(0) / Nx), 2.0)+pow(sin(M_PI * k.e(1) / Ny), 2.0)+pow(sin(M_PI * k.e(2) / Nz), 2.0); 
            kSqu *= pow(2.0 / dx, 2.0);

            if (kSqu > 0.0) {
                std = sqrt(0.5 * constant *
                           exp(-0.5 * kSqu * lc * lc));

		kA[X].gaussian_random(std);
		
            } else {
	      kA[X]=0;
	    }
    }

        FFT_field(kA, A, fft_direction::back);
	
	onsites (ALL)
	  {
	    // if (A[X].norm()>0.0){A[X]=gap*A[X]/A[X].norm();}
	    A[X]=A[X]/sqrt((A[X]*A[X].dagger()).trace());
	  }
	
        pi[ALL] = 0;

        hila::out0 << "k space generation \n";

        break;
  }

  case 2: {
    pi = 0.;                            
    onsites(ALL) {                     
      A[X] = sqrt(0.1) * hila::gaussrand();
    }

    hila::out0 << " normal-phase-real-1 created \n";
    break;    

  }
  case 3: {
    pi = 0.;
    onsites(ALL) {
      /*foralldir(al) foralldir(i){
        A[X].e(al,i).real().gaussian_random();
	}*/
      A[X].gaussian_random();
      A[X] = A[X].real();
    }
    hila::out0 << " normal-phase-real-2 created \n";
    break;    
  }

  case 4: {
    pi = 0.;
    onsites(ALL) {
      foralldir(al) foralldir(i){
	A[X].e(al,i) = sqrt(config.variance_sigma) * hila::gaussian_random<Complex<real_t>>();
      } // doralldir end here
    } // onsites(ALL) end here
  
    hila::out0 << " normal-phase-complex created \n";
    break;
      
  }
    
  case 5: {
    pi = 0;
    real_t gap = MP.gap_B_td(config.Inip, config.IniT);
    hila::out0 <<"Gap B: "<<gap<<"\n";
    onsites(ALL) {
      foralldir(d1)foralldir(d2){

	if (d1==d2){
	  A[X].e(d1,d2).re = 1.0; //hila::gaussrand() hila::random()
	  A[X].e(d1,d2).im = 0.0;
	}
	else {
	  A[X].e(d1,d2).re = 0.0;
	  A[X].e(d1,d2).im = 0.0;}
      }
      A[X] = gap * A[X]/A[X].norm();
    }

    hila::out0 << "Pure B phase \n";

    break;
    }
    
    case 6: {
    pi = 0;
    real_t gapA = MP.gap_A_td(config.Inip, config.IniT);
    real_t gapB = MP.gap_B_td(config.Inip, config.IniT);
    real_t tb = config.Inilc;//config.IniT/ MP.Tcp_mK(config.Inip);
    hila::out0<<"Gap A: "<<gapA<<"\n";
    hila::out0<<"Gap B: "<<gapB<<"\n";
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

    hila::out0 << "Aphase_partial is created \n";

    break;
    }

  case 7: {
    // pi = 0.;
    // real_t gap = MP.gap_A_td(Tp[1], Tp[0]);
    // hila::out0<<"Gap A: "<<gap<<"\n";
    
    // onsites(ALL) {
    //   foralldir(al) foralldir(i){
    // 	A[X].e(al,i) = sqrt(config.variance_sigma) * hila::gaussian_random<Complex<real_t>>();
	
    // 	if ((al==0) && (i==0)) {
    // 	  A[X].e(al,i).re=A[X].e(al,i).re + 1.;
    // 	}
    // 	else if ((al==0) && (i==1)) {
    // 	  A[X].e(al,i).im=A[X].e(al,i).im + 1.;
    // 	} // put bulk A-phase elements into random matrix

    // 	A[X].e(al,i)=(A[X].e(al,i)/sqrt(2.)) * gap;	
    //   } // doralldir end here
    // } // onsites(ALL) end here

    // //A[ALL]=A[x].asArray()
    // hila::out0 << "Aphase_full is created \n";

    break;
  } // case 7: Aphase_full

  case 8: {
    // pi = 0;
    // real_t gapa = MP.gap_A_td(Tp[1], Tp[0]);
    // real_t gapb = MP.gap_B_td(Tp[1], Tp[0]);    
    // // hila::out0<<"Gap A: "<<gap<<"\n";
    // // if (X.coordinate(e_x) == 0 or X.coordinate(e_x) == 1)
    
    // onsites (ALL) {
    // if (
    // 	(X.coordinate(e_x) <= 10 or X.coordinate(e_x) >= (config.lx - 10))
    // 	&& ((X.coordinate(e_y) <= 74 and X.coordinate(e_y) >= 54))
    // 	&& ((X.coordinate(e_z) <= 74 and X.coordinate(e_z) >= 54))
    //    )
    //   {	
    // 	    foralldir(d1)foralldir(d2){
    // 	      if (d1==d2){
    // 		A[X].e(d1,d2).re = 1.0;
    // 		A[X].e(d1,d2).im = 0.0;
    // 	      }
    // 	      else {
    // 		A[X].e(d1,d2).re = 0.0;
    // 		A[X].e(d1,d2).im = 0.0;}
    //           }
    // 	    A[X] = gapb * A[X]/sqrt(3.0);
    //    }
    // //else if (X.coordinate(e_z) == (config.lz - 1) or X.coordinate(e_z) == (config.lz - 2))
    // else if (
    // 	     (X.coordinate(e_x) > 10 and X.coordinate(e_x) < (config.lx - 10))
    // 	     || (
    // 	         (X.coordinate(e_x) <= 10 or X.coordinate(e_x) >= (config.lx - 10))
    // 		 && (!(
    //                     (X.coordinate(e_y) <= 74 and X.coordinate(e_y) >= 54)
    // 	                && (X.coordinate(e_z) <= 74 and X.coordinate(e_z) >= 54)
    //                   ))
    //             )
    // 	    )    
    // 	  {
    // 	    foralldir(d1)foralldir(d2){
    // 	      if (d1==0 && d2==0){
    // 		A[X].e(d1,d2).re = 1.0;
    // 		A[X].e(d1,d2).im = 0.0;
    // 	      }
    // 	      else if (d1==0 && d2==1){
    // 		A[X].e(d1,d2).re = 0.0;  
    // 		A[X].e(d1,d2).im = 1.0;
    // 	      }
    // 	      else {
    // 		A[X].e(d1,d2).re = 0.0;
    // 		A[X].e(d1,d2).im = 0.0;
    // 	      }
    // 	    }
    // 	    A[X] = gapa * A[X]/sqrt(2.0);
    // 	  }
    // }

    // hila::out0 << "B domain in supercooling A \n";

    break;
  } // case 8 block end here
    
  default: {

    // #pragma hila ast_dump
    pi = 0.0; //set derivative matrix to zero
    onsites (ALL) {
      A[X].fill(1.0);
    }
    
    hila::out0 << "Field matrix set to 1 everywhere \n";

    break;
  }
  }


} // initialize() call end here


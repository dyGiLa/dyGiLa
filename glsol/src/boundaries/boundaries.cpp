#define USE_MPI
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>

#include "plumbing/hila.h"
#include "plumbing/fft.h"

#include "glsol.hpp"
#include "matep.hpp"

const phi_t glsol::periodic(phi_t phi1) {

  return phi1;
}

const phi_t glsol::maximal() {

  phi_t phibc;
 
  phibc=0.0;

  return phibc;

}

const phi_t glsol::BphaseBoundary() {


  phi_t phibc;

  real_t gapb = MP.gap_B_td(config.p_boundary, config.T_boundary);
  
  foralldir(d1)foralldir(d2){                                                                                                                                                                     
    if (d1==d2){                                                                                                                                                                                  
      phibc.e(d1,d2).re = 1.0;                                                                                                                                                                     
      phibc.e(d1,d2).im = 0.0;                                                                                                                                                                     
    }                                                                                                                                                                                             
    else {                                                                                                                                                                                        
      phibc.e(d1,d2).re = 0.0;                                                                                                                                                                     
      phibc.e(d1,d2).im = 0.0;}                                                                                                                                                                    
  }                                                                                                                                                                                               

  phibc = gapb * phibc/sqrt(3.0); 
  

  return phibc;

}


const phi_t glsol::AphaseBoundary() {

  phi_t phibc;
  real_t gapa = MP.gap_A_td(config.p_boundary, config.T_boundary);

  foralldir(d1)foralldir(d2){                                                                                                                                                                     
    if (d1==2 && d2==0){                                                                                                                                                                          
      phibc.e(d1,d2).re = 1.0;                                                                                                                                                                     
      phibc.e(d1,d2).im = 0.0;                                                                                                                                                                     
    }                                                                                                                                                                                             
    else if (d1==2 && d2==1){                                                                                                                                                                     
      phibc.e(d1,d2).re = 0.0;  // this A-order parameter same with GL-theory note eq.46                                                                                                           
      phibc.e(d1,d2).im = 1.0;                                                                                                                                                                     
    }                                                                                                                                                                                             
    else {                                                                                                                                                                                        
      phibc.e(d1,d2).re = 0.0;                                                                                                                                                                     
      phibc.e(d1,d2).im = 0.0;                                                                                                                                                                     
    }                                                                                                                                                                                             
  }                                                                                                                                                                                               

  phibc = gapa * phibc/sqrt(2.0);  

  return phibc;
  
}

const phi_t glsol::normalBoundary(){

  
  phi_t phibc;

  foralldir(al) foralldir(i){
      phibc.e(al,i) = sqrt(1.0) * hila::gaussian_random<Complex<real_t>>();
  } // doralldir end here                                                                                                                                                                               
  
  return phibc;
  
}

const phi_t glsol::RobinBoundary(phi_t phiI, int oorN,real_t bt, int face){

  phi_t phibc;
  real_t pbt=config.dx/bt;

  foralldir(d1){
    if(d1==face){
      foralldir(d2){
	phibc.e(d2,d1)=0.0;}
    }else{
      if(oorN == 0){
	foralldir(d2){
	  phibc.e(d2,d1)=phiI.e(d2,d1)*((1.0-pbt)/(1.0+pbt));}
      }else{
	foralldir(d2){
	  phibc.e(d2,d1)=phiI.e(d2,d1)*((1.0+pbt)/(1.0-pbt));}
      }
    }
  }

  return phibc;

}

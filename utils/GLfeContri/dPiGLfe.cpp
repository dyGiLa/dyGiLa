#define USE_MPI 
//#include <sstream>
//#include <iostream>
//#include <iomanip>
//#include <fstream>
//#include <string>
//#include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/fft.h"

#include "glsol.hpp"
#include "matep.hpp"


void glsol::dPiGLfe() {

  // Field<phi_t> deltaPi;
  // Field<Vector<3,Complex<real_t>>> djAaj;  

  /******************************************************/
  /*          Bulk Free energy contribution             */
  /******************************************************/  
  
  onsites (ALL) {
    
    deltaPi[X] = 0; // refresh \delta \Pi

    matep::Matep MPonsites;
    auto AxAt = A[X]*A[X].transpose();
    auto AxAd = A[X]*A[X].dagger();

    real_t beta0 = MPonsites.alpha_td(config.Inip, T[X]);
    real_t beta1 = MPonsites.beta1_td(config.Inip, T[X]);
    real_t beta2 = MPonsites.beta2_td(config.Inip, T[X]);
    real_t beta3 = MPonsites.beta3_td(config.Inip, T[X]);
    real_t beta4 = MPonsites.beta4_td(config.Inip, T[X]);
    real_t beta5 = MPonsites.beta5_td(config.Inip, T[X]);
  
    deltaPi[X] = - beta0*A[X]
      - 2.0*beta1*A[X].conj()*AxAt.trace()
      - 2.0*beta2*A[X]*AxAd.trace()
      - 2.0*beta3*AxAt*A[X].conj()
      - 2.0*beta4*AxAd*A[X]
      - 2.0*beta5*A[X].conj()*A[X].transpose()*A[X]
      - MPonsites.gz_td(config.Inip)*H[X]*(H[X].transpose()*A[X]);

  } // Bulk GL free energy block
  

  /******************************************************/
  /*          Gradient energy contribution              */
  /******************************************************/  

  onsites(ALL) {
    djAaj[X] = 0; // refresh div.A_al    
    foralldir(j) { djAaj[X] += (1.f/(2.*config.dx)) * (A[X + j].column(j) - A[X - j].column(j)); }
  } // div.A_al Block

  onsites(ALL) {
    phi_t didjAaj_X = 0;
    foralldir(di) {
      Vector<3,Complex<real_t>> didjAajVec = djAaj[X + di] - djAaj[X - di];
      for (int al=0; al<NDIM; ++al) didjAaj_X.e(al,di) = didjAajVec[al];
    } // didjAalj matrix block

    deltaPi[X] += (1.f/config.dx) * didjAaj_X;
  } // 2.0 \partial_i djAaj

  onsites (ALL) {
    
    deltaPi[X] += (1.f/(config.dx*config.dx)) * (A[X + e_x] + A[X - e_x]
                                                 + A[X + e_y] + A[X - e_y]
                                                 + A[X + e_z] + A[X - e_z]
                                                 - 6.0*A[X]);
  } // Laplacian of Aalj block

} // dPiGLfe() ends here


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


void glsol::dPiGLfe_AdGRz() {

  const real_t abt_ratio=config.dx/config.bt;

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
  
  onsites(ALL) {
    /******************************************************************/
    /*    start: cook up the A_X+j & A_X-j for AdGRz treatment        */
    /******************************************************************/

    phi_t A_Xmj, A_Xpj;
    AdGRzTreat(A_Xmj, A_Xpj, abt_ratio);
    // if ( X.coordinate(e_z) == 0 )
    //   {
    //     const real_t trcoef=(1.-abt_ratio)/(1.+abt_ratio);
    // 	foralldir(col)
    // 	  {
    //        if(col == 2)
    // 	     { foralldir(row){ A_Xmj.e(row, col)=0.0; } }
    // 	   else
    // 	     { foralldir(row){ A_Xmj.e(row, col)=A[X+e_z].e(row, col)*trcoef;} }
    // 	  } // col loop ends here
    //   }
    // else if ( X.coordinate(e_z) == (config.lz-1) )
    //   {
    //     const real_t trcoef=(1.+abt_ratio)/(1.-abt_ratio);
    // 	foralldir(col)
    // 	  {
    //        if(col == 2)
    // 	     { foralldir(row){ A_Xpj.e(row, col)=0.0; } }
    // 	   else
    // 	     { foralldir(row){ A_Xpj.e(row, col)=A[X-e_z].e(row, col)*trcoef;} }
    // 	  } // col loop ends here
    //   }
    /****************************************************************/
    /*      end: cook up the A_X+j & A_X-j for AdGRz treatment      */
    /****************************************************************/
    
    djAaj[X] = 0;    
    foralldir(j) {
      if ( X.coordinate(e_z) == 0 && j == e_z )
	{ djAaj[X] += (1.f/(2.*config.dx)) * (A[X + j].column(j) - A_Xmj.column(j)); }
      else if ( X.coordinate(e_z) == (config.lz-1) && j == e_z)
	{ djAaj[X] += (1.f/(2.*config.dx)) * (A_Xpj.column(j) - A[X - j].column(j)); }
      else
	{ djAaj[X] += (1.f/(2.*config.dx)) * (A[X + j].column(j) - A[X - j].column(j)); }      
    } // computing div.A_al at X
  } // div.A_al Block, djAalj onsite(ALL) block done

  onsites(ALL) {
    phi_t didjAaj_X = 0;
    foralldir(di) {
      Vector<3,Complex<real_t>> didjAajVec = djAaj[X + di] - djAaj[X - di];
      for (int al=0; al<NDIM; ++al) didjAaj_X.e(al,di) = didjAajVec[al];
    }
    deltaPi[X] += (1.f/config.dx) * didjAaj_X;
  } // 2.0 \partial_i djAaj, djAalj vector field differential block
  
  onsites (ALL) {
    /******************************************************************/
    /* start: cook up the A_X+e_z & A_X-e_z for AdGRz treatment       */
    /******************************************************************/
    phi_t A_Xmez, A_Xpez;
    AdGRzTreat(A_Xmez, A_Xpez, abt_ratio);
    // if ( X.coordinate(e_z) == 0. )
    //   {
    //     const real_t trcoef=(1.-abt_ratio)/(1.+abt_ratio);
    // 	foralldir(col)
    // 	  {
    //        if(col == 2)
    // 	     { foralldir(row){ A_Xmez.e(row, col)=0.0; } }
    // 	   else
    // 	     { foralldir(row){ A_Xmez.e(row, col)=A[X+e_z].e(row, col)*trcoef;} }
    // 	  } // col loop ends here
    //   }
    // else if ( X.coordinate(e_z) == (config.lz-1) )
    //   {
    //     const real_t trcoef=(1.+abt_ratio)/(1.-abt_ratio);
    // 	foralldir(col)
    // 	  {
    //        if(col == 2)
    // 	     { foralldir(row){ A_Xpez.e(row, col)=0.0; } }
    // 	   else
    // 	     { foralldir(row){ A_Xpez.e(row, col)=A[X-e_z].e(row, col)*trcoef;} }
    // 	  } // col loop ends here
    //   }
    /****************************************************************/
    /*   end: cook up the A_X+e_z & A_X-e_z for AdGRz treatment     */
    /****************************************************************/

    if ( X.coordinate(e_z) == 0 )
	{
         deltaPi[X] += (1.f/(config.dx*config.dx))
	               * (A[X + e_x] + A[X-e_x]
                          + A[X + e_y] + A[X-e_y]
                          + A[X + e_z] + A_Xmez
                          - 6.0*A[X]);

	} // starting surface AdGR treatment
    else if ( X.coordinate(e_z) == (config.lz-1) )
	{
         deltaPi[X] += (1.f/(config.dx*config.dx))
	               * (A[X+e_x] + A[X - e_x]
                          + A[X+e_y] + A[X - e_y]
                          + A_Xpez + A[X - e_z]
                          - 6.0*A[X]);
	  
	} // ending surface AdGR treatment
    else
	{
         deltaPi[X] += (1.f/(config.dx*config.dx)) * (A[X + e_x] + A[X - e_x]
                                                      + A[X + e_y] + A[X - e_y]
                                                      + A[X + e_z] + A[X - e_z]
                                                      - 6.0*A[X]);

	}           
  } // Laplacian block done

  /***************************************************************************/
  /*    \delta \pi canonical momentum computing blocks ends from here        */
  /***************************************************************************/
   
} // dPiGLfe_AdGRz() ends here


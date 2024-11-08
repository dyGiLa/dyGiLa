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


void glsol::next_bath() {

  static hila::timer next_timer("timestep");
  Field<phi_t> deltaPi;
  Field<Vector<3,Complex<real_t>>> djAaj;

  Complex<real_t> ep2 = 1.0-exp(-2.0*config.gamma*config.dt);
  real_t tb =  config.IniT/MP.Tcp_mK(config.Inip); // note that the bath tenperature could change. get corretc number

  
  //int bc=config.boundaryConditions;

  next_timer.start();


#ifndef T_FIELD
    real_t gapa = MP.gap_A_td(p, T);
    real_t gapb = MP.gap_B_td(p, T);
#endif
  
  onsites(ALL) {

    A[X] += config.dt * pi[X];

  }

#ifndef T_FIELD
    real_t beta[6];
    point_params(T, p,beta);
#endif
  
  onsites (ALL) {

#ifdef T_FIELD
    real_t beta[6];
    point_params(T[X], p,beta);
#endif
    
    auto AxAt = A[X]*A[X].transpose();
    auto AxAd = A[X]*A[X].dagger();
    phi_t HHA = 0;
    phi_t He = 0;
    real_t sgn=0.0;

    
    deltaPi[X] = - beta[0]*A[X]
      - 2.0*beta[1]*A[X].conj()*AxAt.trace()
      - 2.0*beta[2]*A[X]*AxAd.trace()
      - 2.0*beta[3]*AxAt*A[X].conj()
      - 2.0*beta[4]*AxAd*A[X]
      - 2.0*beta[5]*A[X].conj()*A[X].transpose()*A[X];

  

    foralldir(al) foralldir(i){
      foralldir(k){
	HHA.e(al,i) += H[al]*H[k]*A[X].e(k,i);
      }
    }

    foralldir(al) foralldir(i){
      if (al==i){
	He.e(al,i)=0.0;
      }else {
	if (al<i){
	  sgn=1.0;
	}else{
	  sgn=-1.0;}
	foralldir(k){
	  if (k==al or k==i){
	    He.e(al,i) += 0.0;
	  }else{
	    He.e(al,i) += sgn*H[k];
	  }
	}
      }
    }


    deltaPi[X] -= gh*HHA + gz*He.transpose()*A[X];
    
  }
    
  onsites(ALL) {
    djAaj[X] = 0;
    foralldir(j) {
      //djAaj[X] += A[X + j].column(j) - A[X - j].column(j);
      djAaj[X] += shift(A[X+j],A[X],A[X-j], X.coordinates(),j,1).column(j)-shift(A[X+j],A[X],A[X-j], X.coordinates(),j,0).column(j);
    }
  }

  onsites(ALL) {
    phi_t mat;
    foralldir(d) {
      auto col = djAaj[X+d] - djAaj[X-d];
      for (int i=0; i<NDIM; i++) mat.e(i,d) = col[i];
    }

    deltaPi[X] += (1.0/(2.0*(config.dx*config.dx)))*mat;
  }

  onsites (ALL) {

    deltaPi[X] +=  (1.0 / (config.dx * config.dx)) * (shift(A[X+e_x],A[X], A[X-e_x],X.coordinates(),e_x,1) + shift(A[X+e_x],A[X], A[X-e_x],X.coordinates(),e_x,0)//A[X + e_x] + A[X - e_x]
						      + shift(A[X+e_y],A[X], A[X-e_y],X.coordinates(),e_y,1) + shift(A[X+e_y],A[X], A[X-e_y],X.coordinates(),e_y,0) //+ A[X + e_y] + A[X - e_y]
						      + shift(A[X+e_z],A[X], A[X-e_z],X.coordinates(),e_z,1) + shift(A[X+e_z],A[X], A[X-e_z],X.coordinates(),e_z,0) //+ A[X + e_z] + A[X - e_z]
						      - 6.0 * A[X]);
    
    //deltaPi[X] += (1.0/(4.0*config.dx*config.dx)) * (A[X + e_x] + A[X - e_x]
    //                                               + A[X + e_y] + A[X - e_y]
    //                                               + A[X + e_z] + A[X - e_z]
    //                                               - 6.0*A[X]);

  }

  if (t < config.tdif)
    {
      pi[ALL] = deltaPi[X]/(config.difFac);
      t += config.dt/config.difFac;
    }
  else if (t < config.tdis && config.gamma.squarenorm() > 0 )
    {
      onsites(ALL){
	phi_t rad_mat;
	rad_mat.gaussian_random();
	if ((config.useTbath == 1)
	    && (t >= config.Tbath_start))
	  {
	    pi[X] = pi[X] + (deltaPi[X] - 2.0 * config.gamma * pi[X])*(config.dt/2.0);
	    pi[X] = sqrt(1.0-ep2)*pi[X] + sqrt(ep2)*tb*rad_mat;
	    pi[X] = pi[X] + (deltaPi[X] - 2.0 * config.gamma * pi[X])*(config.dt/2.0);
	  }
	else
	  {
	    pi[X] = pi[X] + (deltaPi[X] - 2.0 * config.gamma * pi[X])*(config.dt/2.0);
	  }
      }
      
      t += config.dt;

    }
  else
    {
      pi[ALL] = pi[X] + deltaPi[X]*config.dt;
      t += config.dt;
    }

  //hila::out0<<"Per mod: "<<modP / lattice.volume()<<"/n";
  next_timer.stop();

} // next_bath() ends here


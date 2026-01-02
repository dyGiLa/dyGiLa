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
  // Field<phi_t> deltaPi;
  // Field<Vector<3,Complex<real_t>>> djAaj;

  const real_t Tcp_mK = MP.Tcp_mK(config.Inip);
  // const real_t kBTCf0p_ratio = MP.kBTCf0p_ratio(config.Inip);
  // const real_t volElemLattice = config.dx * config.dx * config.dx;
  
  int bc=config.boundaryConditions;
  // hila::out0 << "bc is " << bc << " in this next_bath() call " << std::endl;

  next_timer.start();

  onsites(ALL) {

    matep::Matep MPonsites;
    
    real_t gapa = MPonsites.gap_A_td(config.Inip, T[X]);
    real_t gapb = MPonsites.gap_B_td(config.Inip, T[X]);

    A[X] += config.dt * pi[X];

    if (bc == 1)
      {
        if (X.coordinate(e_z) == 0 or X.coordinate(e_z) == 1)
          {
            foralldir(d1)foralldir(d2){
              if (d1==d2){
                A[X].e(d1,d2).re = 1.0;
                A[X].e(d1,d2).im = 0.0;
              }
              else {
                A[X].e(d1,d2).re = 0.0;
                A[X].e(d1,d2).im = 0.0;
	      }	
            }
            A[X] = gapb * A[X]/sqrt(3.0);
          }
        else if (X.coordinate(e_z) == (config.lz - 1) or X.coordinate(e_z) == (config.lz - 2))
          {
            foralldir(d1)foralldir(d2){
              if (d1==2 && d2==0){
                A[X].e(d1,d2).re = 1.0;
                A[X].e(d1,d2).im = 0.0;
              }
              else if (d1==2 && d2==1){
                A[X].e(d1,d2).re = 0.0;
                A[X].e(d1,d2).im = 1.0;
              }
              else {
                A[X].e(d1,d2).re = 0.0;
                A[X].e(d1,d2).im = 0.0;
              }
	               }
            A[X] = gapa * A[X]/sqrt(2.0);
          }
        }
    else if (bc == 2)
      {
        if (
	    X.coordinate(e_x) == 0 || X.coordinate(e_x) == (config.lx - 1) ||
            X.coordinate(e_x) == 1 || X.coordinate(e_x) == (config.lx - 2) ||
            X.coordinate(e_y) == 0 || X.coordinate(e_y) == (config.ly - 1) ||
            X.coordinate(e_y) == 1 || X.coordinate(e_y) == (config.ly - 2) ||
            X.coordinate(e_z) == 0 || X.coordinate(e_z) == (config.lz - 1) ||
            X.coordinate(e_z) == 1 || X.coordinate(e_z) == (config.lz - 2)
	   )
          {
            A[X]=0.0;	    
          }
      }
  } // onsite() block ends here

  dPiGLfe(); // compute free energy contribution for delta Pi
  ABOBA();   // Canonical momentum Langevin update
  
  // if (t < config.tdif)
  //   {
  //     pi[ALL] = deltaPi[X]/(config.difFac);
  //     t += config.dt/config.difFac;
  //   }
  // else if (t < config.tdis && config.useTbath == 1)
  //   {
  //     onsites(ALL){
  //       matep::Matep MPonsites;

  //       real_t ep2 = 1.0-exp(-2.0 * MPonsites.gamma_td(config.Inip, T[X], phaseMarker[X]) * config.dt);

  // 	phi_t rad_mat;
  // 	rad_mat.gaussian_random();
	
  //       // damping term gives 2.0, but it is absobed by new defination of gamma, then coef is 1.0	
  // 	pi[X] = pi[X] + (deltaPi[X] - 1.0 * MPonsites.gamma_td(config.Inip, T[X], phaseMarker[X]) * pi[X])*(config.dt/2.0);
	
  //       pi[X] = sqrt(1.0 - ep2) * pi[X] + sqrt(ep2 * (T[X]/Tcp_mK) * (kBTCf0p_ratio/volElemLattice)) * rad_mat;

  //       // damping term gives 2.0, but it is absobed by new defination of gamma, then coef is 1.0      
  //       pi[X] = pi[X] + (deltaPi[X] - 1.0 * MPonsites.gamma_td(config.Inip, T[X], phaseMarker[X]) * pi[X])*(config.dt/2.0);
	
  //     }

  //     t += config.dt;
  //   }
  // else
  //   {
  //     pi[ALL] = pi[X] + deltaPi[X]*config.dt;
  //     t += config.dt;
  //   }

  //hila::out0<<"Per mod: "<<modP / lattice.volume()<<"/n";
  next_timer.stop();

} // next_bath() ends here


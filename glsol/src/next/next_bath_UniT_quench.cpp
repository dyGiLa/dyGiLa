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


void glsol::next_bath_UniT_quench() {

  static hila::timer next_timer("timestep");
  // Field<phi_t> deltaPi;
  // Field<Vector<3,Complex<real_t>>> djAaj;

  const real_t Tcp_mK_Inip = MP.Tcp_mK(config.Inip);
  // const real_t kBTCf0p_ratio = MP.kBTCf0p_ratio(config.Inip);
  // const real_t volElemLattice = config.dx * config.dx * config.dx;
  
  int bc=config.boundaryConditions;
  // hila::out0 << "bc is " << bc << " in this next_bath() call " << std::endl;

  next_timer.start();

  std::initializer_list<int> coordsList {0,0,0};
  const CoordinateVector originpoints(coordsList);
  if ( T.get_element(originpoints) > (config.Ttd_Qend * MP.Tcp_mK(config.Inip)) )
    {
     if (
	 (config.has1stQStop == true)
	 && ((T.get_element(originpoints) - (config.Ttd_Q1st * MP.Tcp_mK(config.Inip))) <= 0.0)
	 && (t < config.tQ1Waiting)
	)
       {/*empty block*/}
     else
       {
	// 1st quench and 2nd quench has different tauQ
	config.tauQ = (t > config.tQ1Waiting) ? config.tauQ2 : config.tauQ1;
	 
        // Temperature update for uniform quench
        T[ALL] = T[X] - ((config.dt/config.tauQ) * Tcp_mK_Inip);
        // hila::out0 << " T in site is " << T.get_element(originpoints) << std::endl;

       }
    }
      
  onsites(ALL) {

    matep::Matep MPonsites;
    real_t gapa = MPonsites.gap_A_td(config.Inip, T[X]);
    real_t gapb = MPonsites.gap_B_td(config.Inip, T[X]);

    A[X] += config.dt * pi[X];

    if (bc == 1) // you don't want use this, especially on GPU
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
  ABOBA();   // canonicl momentum Langevin update
  // if (t < config.tdif)
  //   {
  //     pi[ALL] = deltaPi[X]/(config.difFac);
  //     t += config.dt/config.difFac;
  //   }
  // else if (t < config.tdis && config.gamma.squarenorm() > 0 )
  //   {

  //     onsites(ALL){
	
  // 	phi_t rad_mat;       	
  // 	rad_mat.gaussian_random();

  //       matep::Matep MPonsites;

  //       real_t ep2 = 1.0-exp(-2.0 * MPonsites.gamma_td(config.Inip, T[X], phaseMarker[X]) * config.dt);	

  // 	// damping term gives 2.0, but it is absobed by new defination of gamma, then coef is 1.0	
  // 	pi[X] = pi[X] + (deltaPi[X] - 1.0 * MPonsites.gamma_td(config.Inip, T[X], phaseMarker[X]) * pi[X])*(config.dt/2.0);

  // 	//pi[X] = sqrt(1.0-ep2)*pi[X] + sqrt(ep2)*tb*rad_mat;
  // 	pi[X] = sqrt(1.0-ep2)*pi[X] + sqrt(ep2 * (T[X]/Tcp_mK_Inip) * (kBTCf0p_ratio/volElemLattice))*rad_mat;

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


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


void glsol::next_UniT_Hfield_constrained() {

  static hila::timer next_timer("timestep");
  Field<phi_t> deltaPi;
  Field<Vector<3,Complex<real_t>>> djAaj;

  //real_t tb =  config.IniT/MP.Tcp_mK(config.Inip);

  //hila::out0 <<"Bath evolution with: ep2="<<ep2<<" and tb="<<tb<<"\n";

  //double modP=0.0;
  
  int bc=config.boundaryConditions;
  hila::out0 << "bc is " << bc << " in this next_bath() call " << std::endl;

  next_timer.start();

  std::initializer_list<int> coordsList {0,0,0};
  const CoordinateVector originpoints(coordsList); 
      
  onsites(ALL) {

    real_t gap_A = MP.gap_A_td(config.Inip, T[X]);
    real_t gap_B = MP.gap_B_td(config.Inip, T[X]);

    A[X] += config.dt * pi[X];

    if (bc == 1)
      {
        if (
	    (X.coordinate(e_x) == 0)
	    || (X.coordinate(e_x) == 1)
	   )
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
            A[X] = gap_B * A[X]/sqrt(3.0);
          }
        else if (
		 (X.coordinate(e_x) == (config.lx - 1))
		 || (X.coordinate(e_x) == (config.lx - 2))
		)
          {
            foralldir(d1)foralldir(d2){
              if (d1==0 && d2==0){
                A[X].e(d1,d2).re = 1.0;
                A[X].e(d1,d2).im = 0.0;
              }
              else if (d1==0 && d2==1){
                A[X].e(d1,d2).re = 0.0;
                A[X].e(d1,d2).im = 1.0;
              }
              else {
                A[X].e(d1,d2).re = 0.0;
                A[X].e(d1,d2).im = 0.0;
              }
	    }
            A[X] = gap_A * A[X]/sqrt(2.0);
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

  // evaluate the volume avraged int_v A^dagger A
  const double V = lattice.volume();
  Complex<double> sumtrAdA{0.}, ViinttrAdA{0.};

  //hila::set_allreduce(false);
  onsites(ALL) {
    sumtrAdA += ((A[X].dagger())*A[X]).trace();
  }

  const double gap_A = MP.gap_A_td(config.Inip, config.IniT),
               gap_B = MP.gap_B_td(config.Inip, config.IniT);
  
  ViinttrAdA  = sumtrAdA/V;
  double constrained_coef = ViinttrAdA.re
                            - ((1.-config.kappa)*gap_A*gap_A + config.kappa*gap_B*gap_B); 
  hila::out0 << " V^-1 int_V intAdA is " << ViinttrAdA
             << ", constrained_coef is " << constrained_coef
             << "\n"
	     << std::endl;

  const real_t lambda = (t > config.confSmoothTime)
                        ? config.lambda1
                        : config.lambda0;
  
  onsites (ALL) {

    real_t beta[6];
    point_params(T[X], config.Inip, beta);

    auto AxAt = A[X]*A[X].transpose();
    auto AxAd = A[X]*A[X].dagger();

    deltaPi[X] = - beta[0]*A[X]
      - 2.0*beta[1]*A[X].conj()*AxAt.trace()
      - 2.0*beta[2]*A[X]*AxAd.trace()
      - 2.0*beta[3]*AxAt*A[X].conj()
      - 2.0*beta[4]*AxAd*A[X]
      - 2.0*beta[5]*A[X].conj()*A[X].transpose()*A[X]
      - MP.gz_td(config.Inip)*H[X]*(H[X].transpose()*A[X])
      /* the following term is constrained term contributions */      
      - 2.0*lambda*constrained_coef*A[X];
  } 

  onsites(ALL) {
    djAaj[X] = 0;
    foralldir(j) {
      djAaj[X] += A[X + j].column(j) - A[X - j].column(j);
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
    
    deltaPi[X] += (1.0/(4.0*config.dx*config.dx)) * (A[X + e_x] + A[X - e_x]
                                                     + A[X + e_y] + A[X - e_y]
                                                     + A[X + e_z] + A[X - e_z]
                                                     - 6.0*A[X]);
  }

    //onsites (ALL) {deltaPi[X] *= config.dt;} // I think that this is the problem, multiplication with respect to dt   
  if (t < config.tdif)
    {
      pi[ALL] = deltaPi[X]/(config.difFac);
      t += config.dt/config.difFac;
    }
  else if (t < config.tdis && config.gamma.squarenorm() > 0 )
    {
      pi[ALL] = pi[X] + (deltaPi[X] - 2.0 * config.gamma * pi[X])*(config.dt);

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


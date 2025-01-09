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


void glsol::next_bath_hotblob_quench_Hfield() {

  static hila::timer next_timer("timestep");
  Field<phi_t> deltaPi;
  Field<Vector<3,Complex<real_t>>> djAaj;

  //real_t tb =  config.IniT/MP.Tcp_mK(config.Inip);

  //hila::out0 <<"Bath evolution with: ep2="<<ep2<<" and tb="<<tb<<"\n";

  //double modP=0.0;
  
  int bc=config.boundaryConditions;
  // hila::out0 << "bc is " << bc << " in this next_bath() call " << std::endl;

  next_timer.start();

  // std::initializer_list<int> coordsList {0,0,0};
  // const CoordinateVector originpoints(coordsList);

  // update the Temperature field
  // compute new blob profile on next time step
  if (t > config.tStats /* tStats should be 0 */)
    {
      onsites(ALL)
	{
	  /* hila's coordinate index is counted from zero at corner,
           * so for a blob at the center of box, you need coordinate transformation
           */
          auto x = X.coordinate(e_x) - config.lx/2.0;
	  auto y = X.coordinate(e_y) - config.ly/2.0;
	  auto z = X.coordinate(e_z) - config.lz/2.0;

          real_t tm  = MP.t_TcMax_blob(config.Inip, config.Ttdb1, config.Ttdb0, config.t1);
	  real_t Tcp_mK = MP.Tcp_mK(config.Inip);

	  auto r2 = x*x + y*y + z*z;

	  T[X] = ((config.Ttdb1 - config.Ttdb0) * Tcp_mK
		  * std::pow((config.t1/(tm + t)), 3./2.)
		  * exp(-r2/(4. * MP.Dd(config.Inip) * (tm + t))))
	         + config.Ttdb0 * Tcp_mK;
	} // onsites(ALL) block ends here
    }
  
  // hila::out0 << " T in site is " << T.get_element(originpoints) << std::endl;

  // update the random weight in Langevin eqn with updated uniform gamma 
  // Complex<real_t> ep2 = 1.0-exp(-2.0*config.gamma*config.dt);
  

  
  onsites(ALL) {

    real_t gapa = MP.gap_A_td(config.Inip, T[X]);
    real_t gapb = MP.gap_B_td(config.Inip, T[X]);

    A[X] += config.dt * pi[X];

    if (bc == 1)
      {/* empty block */}
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
      - MP.gz_td(config.Inip)*H[X]*(H[X].transpose()*A[X]);

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

    //Matep MP;
    //real_t tb = config.IniT/ MP.Tcp_mK(config.Inip);
    //real_t sig = sqrt(2.0*tb*config.gamma); //should we have t
    //phi_t rad_mat;
    
    deltaPi[X] += (1.0/(4.0*config.dx*config.dx)) * (A[X + e_x] + A[X - e_x]
                                                     + A[X + e_y] + A[X - e_y]
                                                     + A[X + e_z] + A[X - e_z]
                                                     - 6.0*A[X]);
    //rad_mat.gaussian_random();
    //deltaPi[X] += rad_mat*sig;

  }

    //onsites (ALL) {deltaPi[X] *= config.dt;} // I think that this is the problem, multiplication with respect to dt   
  if (t < config.tdif)
    {
      pi[ALL] = deltaPi[X]/(config.difFac);
      t += config.dt/config.difFac;
    }
  else if (t < config.tdis && config.useTbath == 1 )
    {
      onsites(ALL){

	//Complex<real_t> ep2 = 1.0-exp(-2.0 * MP.gamma_td(config.Inip, T[X]) * config.dt);
	real_t ep2 = 1.0-exp(-2.0 * MP.gamma_td(config.Inip, T[X], phaseMarker[X]) * config.dt);
	
	phi_t rad_mat;
	rad_mat.gaussian_random();
	
	// damping term gives 2.0, but it is absobed by new defination of gamma, then coef is 1.0	
	pi[X] = pi[X] + (deltaPi[X] - 1.0 * MP.gamma_td(config.Inip, T[X], phaseMarker[X]) * pi[X])*(config.dt/2.0);	

	//pi[X] = sqrt(1.0-ep2)*pi[X] + sqrt(ep2)*tb*rad_mat;
	/* Langevin refresh  */
	pi[X] = sqrt(1.0 - ep2) * pi[X] + sqrt(ep2) * (T[X]/MP.Tcp_mK(config.Inip)) * rad_mat; 
	//modP += sqrt(ep2)*tb*rad_mat.norm();
      }

      // damping term gives 2.0, but it is absobed by new defination of gamma, then coef is 1.0	
      pi[ALL] = pi[X] + (deltaPi[X] - 1.0 * MP.gamma_td(config.Inip, T[X], phaseMarker[X]) * pi[X])*(config.dt/2.0);      

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


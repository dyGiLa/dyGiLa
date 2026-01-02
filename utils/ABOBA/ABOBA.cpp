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


void glsol::ABOBA() {

  /******************************************************/
  /*      Canonical momentum Langevin update            */
  /******************************************************/  

  const real_t Tcp_mK = MP.Tcp_mK(config.Inip);
  const real_t kBTCf0p_ratio = MP.kBTCf0p_ratio(config.Inip);
  const real_t volElemLattice = config.dx * config.dx * config.dx;
  
  if (t < config.tdis && config.useTbath == 1 )
    {
      onsites(ALL){
	matep::Matep MPonsites;

	//Complex<real_t> ep2 = 1.0-exp(-2.0 * MP.gamma_td(config.Inip, T[X]) * config.dt);
	real_t ep2 = 1.0-exp(-2.0 * MPonsites.gamma_td(config.Inip, T[X], phaseMarker[X]) * config.dt);
	
	phi_t rad_mat;
	rad_mat.gaussian_random();
	
	// damping term gives 2.0, but it is absobed by new defination of gamma, then coef is 1.0	
	pi[X] = pi[X] + (deltaPi[X] - 1.0 * MPonsites.gamma_td(config.Inip, T[X], phaseMarker[X]) * pi[X])*(config.dt/2.0);	

	//pi[X] = sqrt(1.0-ep2)*pi[X] + sqrt(ep2)*tb*rad_mat;
	/* Langevin refresh  */
	pi[X] = sqrt(1.0 - ep2) * pi[X] + sqrt(ep2 * (T[X]/Tcp_mK) * (kBTCf0p_ratio/volElemLattice)) * rad_mat; 
	//modP += sqrt(ep2)*tb*rad_mat.norm();

        // damping term gives 2.0, but it is absobed by new defination of gamma, then coef is 1.0	
        pi[X] = pi[X] + (deltaPi[X] - 1.0 * MPonsites.gamma_td(config.Inip, T[X], phaseMarker[X]) * pi[X])*(config.dt/2.0);      
	
      }

      t += config.dt;
    }
  else
    {
      pi[ALL] = pi[X] + deltaPi[X]*config.dt;
      t += config.dt;
    }

} // ABOBA() ends here


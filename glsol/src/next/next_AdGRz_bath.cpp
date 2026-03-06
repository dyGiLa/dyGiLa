#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/fft.h"

#include "glsol.hpp"
#include "matep.hpp"


void glsol::next_AdGRz_bath() {

  static hila::timer next_timer("timestep");
  // Field<phi_t> deltaPi;
  // Field<Vector<3,Complex<real_t>>> djAaj;

  // const real_t Tcp_mK = MP.Tcp_mK(config.Inip);
  // const real_t kBTCf0p_ratio = MP.kBTCf0p_ratio(config.Inip);
  // const real_t volElemLattice = config.dx * config.dx * config.dx;
  
  // int bc=config.boundaryConditions;

  next_timer.start();

  /***************************************************************************/
  /*  block for OP A field update and normal component Dirichlet BC          */
  /***************************************************************************/
  onsites (ALL)
    {
      // bulk update      
      A[X] += config.dt * pi[X];
      
      if (X.coordinate(e_z) == 0 || X.coordinate(e_z) == (config.lz - 1))
          {
	    A[X].e(0,2).re=0.; A[X].e(0,2).im=0.;
            A[X].e(1,2).re=0.; A[X].e(1,2).im=0.;
            A[X].e(2,2).re=0.; A[X].e(2,2).im=0.;	    	    
	  }  

    } // AdGRz Aal_z pair-breaking BC
  
  dPiGLfe_AdGRz();
  ABOBA_gBranch();  
  
  next_timer.stop();

} // next_bath() ends here


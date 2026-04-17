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


void glsol::next_bath_Quasi2Dhotblob_quench_AdGR_Hfield() {

  static hila::timer next_timer("timestep");
  // Field<phi_t> deltaPi;
  // Field<Vector<3,Complex<real_t>>> djAaj;

  const real_t Tcp_mK = MP.Tcp_mK(config.Inip);
  // const real_t kBTCf0p_ratio = MP.kBTCf0p_ratio(config.Inip);
  // const real_t volElemLattice = config.dx * config.dx * config.dx;
  
  // int bc=config.boundaryConditions;
  // hila::out0 << "bc is " << bc << " in this next_bath() call " << std::endl;

  next_timer.start();

  // std::initializer_list<int> coordsList {0,0,0};
  // const CoordinateVector originpoints(coordsList);

  // update the Temperature field
  // compute new blob profile on next time step

  real_t tm = MP.t_TcMax_blob(config.Inip, config.Ttdb1, config.Ttdb0, config.t1);
  real_t t0 = MP.t_TcVanish_blob(config.Inip, config.Ttdb1, config.Ttdb0, config.t1);
  
  if (t > config.tStats /* tStats should be > 0 */)
    {
      onsites(ALL)
	{
	  /* hila's coordinate index is counted from zero at corner,
           * so for a blob at the center of box, you need coordinate transformation
           */
    // auto x = X.coordinate(e_x) - config.lx/2.0;
	  // auto y = X.coordinate(e_y) - config.ly/2.0;
	  // auto z = X.coordinate(e_z) - config.lz/2.0;

	  matep::Matep MPonsites;

	  // auto r2 = (x*x + y*y + z*z)/4.f;
    auto r2 = (X.coordinates() - lattice.size() / 2).squarenorm() * sqr(config.dx);

    if (config.Blob_Tc_cutoff == true)
	    {
	     /* conduct Tc cutoff  */ 
             if ( (t + tm) < t0 )
	       {		 
	        // compute Tc frontier after time tm, you want this stays inside if block
                real_t r2Tc = MPonsites.r2_Tc_blob(config.Inip, config.use_CustomerDctxi, config.Dctxi, config.Ttdb1, config.Ttdb0, config.t1, t + tm);  

	        // t1 is in unit of tGL
	        T[X] = (r2 <= r2Tc)
	               ? ((config.Ttdb1 - config.Ttdb0) * Tcp_mK
	   	          * std::pow((config.t1/(tm + t)), 3./2.)
		          * exp(-r2/(4. * MPonsites.Dd(config.Inip, config.use_CustomerDctxi, config.Dctxi) * (tm + t))))
	                 + config.Ttdb0 * Tcp_mK
	               : config.Ttdb0 * Tcp_mK;
	       }
	     else
	       T[X] = config.Ttdb0 * Tcp_mK;
	    } // Tc cutoff block ends here
	  else
	    /* Normal phase D and C computed A- & Normal-phase T profile */
	    T[X] = ((config.Ttdb1 - config.Ttdb0) * Tcp_mK
	   	    * std::pow((config.t1/(tm + t)), 3./2.)
		    * exp(-r2/(4. * MPonsites.Dd(config.Inip, config.use_CustomerDctxi, config.Dctxi) * (tm + t))))
		   + config.Ttdb0 * Tcp_mK;
	    	 	  	  
	} // onsites(ALL) block ends here
    }

  /***************************************************************************/
  /*         OP A field update & normal component Dirichlet BC               */
  /***************************************************************************/
  
  onsites(ALL) {
    // matep::Matep MPonsites;
    
    // real_t gapa = MPonsites.gap_A_td(config.Inip, T[X]);
    // real_t gapb = MPonsites.gap_B_td(config.Inip, T[X]);

    A[X] += config.dt * pi[X];

    if (X.coordinate(e_z) == 0 || X.coordinate(e_z) == (config.lz - 1))
        {
         A[X].e(0,2).re=0.; A[X].e(0,2).im=0.;
         A[X].e(1,2).re=0.; A[X].e(1,2).im=0.;
         A[X].e(2,2).re=0.; A[X].e(2,2).im=0.;	    	    
	}      
  } // onsite() block ends here

  /***************************************************************************/
  /*      OP A field update & normal component Dirichlet BC ends             */
  /***************************************************************************/

  dPiGLfe_AdGRz();
  ABOBA_gBranch();
    
  next_timer.stop();

} // next_bath() ends here


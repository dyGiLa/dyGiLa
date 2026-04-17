#define USE_MPI 
// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
#include <string>
//#include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/fft.h"

#include "glsol.hpp"
#include "matep.hpp"


void glsol::next_bath_UniT_quench_AdGRz_Hfield() {

  static hila::timer next_timer("timestep");
  // Field<phi_t> deltaPi;
  // Field<Vector<3,Complex<real_t>>> djAaj;

  const real_t Tcp_mK = MP.Tcp_mK(config.Inip);
  // const real_t kBTCf0p_ratio = MP.kBTCf0p_ratio(config.Inip);
  // const real_t volElemLattice = config.dx * config.dx * config.dx;    
  int bc=config.boundaryConditions;

  next_timer.start();

  std::initializer_list<int> coordsList {0,0,0};
  const CoordinateVector originpoints(coordsList);

  // update the Temperature field
  if (
       //------------------------------------------------------------
       // Here one can have quench and anti-quench
       // i.e., cooling or warming
       //------------------------------------------------------------      
       (T.get_element(originpoints) > (config.Ttd_Qend * MP.Tcp_mK(config.Inip))
	&& config.use_antiQuench == false)
       ||
       (T.get_element(originpoints) < (config.Ttd_Qend * MP.Tcp_mK(config.Inip))
	&& config.use_antiQuench == true)       
     )
    {
     if (
	 //------------------------------------------------------------
	 // two stage homogenuous quench block,
	 // set has1stQStop be 0 do the stright quench
	 //------------------------------------------------------------	 	 
	 ((config.has1stQStop == true)
	  && ((T.get_element(originpoints) - (config.Ttd_Q1st * MP.Tcp_mK(config.Inip))) <= 0.0)
	  && (t < config.tQ1Waiting)
	  && (config.use_antiQuench == false))
	 ||
	 ((config.has1stQStop == true)
	  && ((T.get_element(originpoints) - (config.Ttd_Q1st * MP.Tcp_mK(config.Inip))) >= 0.0)
	  && (t < config.tQ1Waiting)
	  && (config.use_antiQuench == true))	 
	)
       {/*empty block*/}
     else
       {
	//------------------------------------------------------------
	// the 2nd homogenuous quench block 
	// 1st quench and 2nd quench may have different tauQ
	// however, if tauQ1 in fact equal to tauQ2, as well as config.has1stQStop == false
	// calling this functionequal to homogenous quench with one tauQ continously
        //------------------------------------------------------------	 	 
	config.tauQ = (t > config.tQ1Waiting) ? config.tauQ2 : config.tauQ1;        	 
        // Temperature update for uniform quench
	onsites(ALL)
	  {
	    T[X] = ( config.use_antiQuench == true ) ? T[X] + ((config.dt/config.tauQ) * Tcp_mK)
	                                             : T[X] - ((config.dt/config.tauQ) * Tcp_mK);
	  }
        // hila::out0 << " T in site is " << T.get_element(originpoints) << std::endl;
       }
    } // Temperature handling block ends here

  /***************************************************************************/
  /*         OP A field update & normal component Dirichlet BC               */
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

  /***************************************************************************/
  /*      OP A field update & normal component Dirichlet BC ends             */
  /***************************************************************************/
  
  dPiGLfe_AdGRz();
  ABOBA_gBranch();
  
  next_timer.stop();

} // next_bath() ends here


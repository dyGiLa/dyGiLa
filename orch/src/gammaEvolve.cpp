//#define USE_PARIO 
#define USE_MPI 
// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <string>
// #include <assert.h>

// #include "plumbing/hila.h"
// #include "plumbing/globals.h" 

#include "glsol.hpp"
// #include "matep_namespace_utils.hpp"
#include "orch.hpp"

// #if defined USE_PARIO 
// #include "pario.hpp"
// #endif

namespace orch {
  
void gammaEvolve(glsol &gl, unsigned int &stat_counter, const unsigned int &steps, const CoordinateVector &originpoints) {
  /* config.gamma update block */
  if (
      (gl.config.TDependnetgamma == 1)
      && (gl.config.initialConditionT != 2)
      && (stat_counter < (gl.config.gammaoffc)*steps)
     )
     // update gl.config.gamma if T-dependency of gamma is turned on
   {
     gl.config.gamma = (gl.T.get_element(originpoints) < gl.MP.Tcp_mK(gl.config.Inip))
                       ? gl.MP.gamma_td(gl.config.Inip, gl.T.get_element(originpoints), gl.phaseMarker.get_element(originpoints))
		       : gl.MP.gamma_td(gl.config.Inip, gl.MP.Tcp_mK(gl.config.Inip), gl.phaseMarker.get_element(originpoints));
		
   }
  else if (stat_counter >= (gl.config.gammaoffc)*steps)
      	  // Set gamma to 2nd value after gamma switching count for all profile,
          // this is for configuraton capture.
          { gl.config.gamma = gl.config.gamma2; }	      	    
  
} /* config.gamma update block end here */

} // orch namespace ends here


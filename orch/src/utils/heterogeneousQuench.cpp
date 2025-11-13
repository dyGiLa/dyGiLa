//#define USE_PARIO 
#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
//#include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/globals.h" 

#include "glsol.hpp"
//#include "matep_namespace_utils.hpp"
#include "orch.hpp"
#include "orch_utils.hpp"

// #if defined USE_PARIO 
// #include "pario.hpp"
// #endif

namespace orch {
  
void heterogeneousQuench(glsol &gl, unsigned int &stat_counter, const unsigned int &steps, const CoordinateVector &originpoints, const unsigned int &modSteps, const unsigned int &modPSSR)
{
    if (stat_counter < (gl.config.gammaoffc)*steps)
      {
       //hila::out0 << "just before call next-blob() " << std::endl;
       if ( gl.config.use_AdGRz_surfaces == 1 )
	 // Quasi 2D Cylinderial blob with AdGR boundary for high-T
	 { /*gl.next_bath_Quasi2Dhotblob_quench_AdGR_Hfield();*/ }
       else
	 // 3D blob without boundary surface effect
         { gl.next_bath_hotblob_quench_Hfield(); }
       
         if (
             (modSteps == 0)
	     && (modPSSR == 0)
            )
           {// squeze IO a little bit
            hila::out0 << "gl.t is " << gl.t 
	               << ", next_bath_hotblob_quench_Hfield() call " 
	               << ", Tc is " << gl.MP.Tcp_mK(gl.config.Inip)
	               << ", Ttdb0 is " << gl.config.Ttdb0	      
                       << std::endl;
           }

      } // blob evolution block
    else if (stat_counter >= (gl.config.gammaoffc)*steps)
      { // configuration catch block
       ++gl.extinguish_t; //estinguish time count, in step of dt
       gl.next_bath_hotblob_quench_Hfield_confCatch();
       if (
           (modSteps == 0)
	   && (modPSSR == 0)	   
          )
         {// squeze IO a little bit
          hila::out0 << "gl.t is " << gl.t 
	             << ", next_bath_hotblob_quench_Hfield_confCatch() call " 
	             << ", Tc is " << gl.MP.Tcp_mK(gl.config.Inip)
	             << ", Ttdb0 is " << gl.config.Ttdb0	      
                     << std::endl;
         }

      } // conf catch block
    
} // heterogenous quench, hot blob T-profile 

} // orch namespace ends here

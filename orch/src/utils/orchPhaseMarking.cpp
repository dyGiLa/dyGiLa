//#define USE_PMD_GAMMA
//#define USE_PARIO 
#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/globals.h" 

#include "glsol.hpp"
//#include "matep_namespace_utils.hpp"
//#include "orch.hpp"
#include "orch_utils.hpp"

// #if defined USE_PARIO 
// #include "pario.hpp"
// #endif

namespace orch {
  
  void phaseMarking(glsol &gl, unsigned int &stat_counter, const unsigned int &steps)
   {
#ifdef USE_PMD_GAMMA
     // do every-dt phase-Marking when gamma is T-dependent heterogenously	  
     gl.phaseMarking();
     if (stat_counter % steps == 0)     
#endif	    
#ifndef USE_PMD_GAMMA	      
       gl.phaseMarking();
#endif	    	      
       if ((stat_counter / steps) % gl.config.PSSRatio == 0)
	 { hila::out0 << "gl.t is " << gl.t
		      << ", phaseMarking() call is done. "
		      << std::endl; }
  
   } // orch::phaseMarking func ends here

} // orch namespace ends here

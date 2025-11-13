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
  
void nextBlocks(glsol &gl, unsigned int &stat_counter, const unsigned int &steps, const CoordinateVector &originpoints)
{
  /*******************************************************************/
  /* the following if else blocks are different system synamic updates  
   *  all next_xxx() functions call happen here.
   */
  /*******************************************************************/

  const unsigned int modSteps = stat_counter % steps;
  const unsigned int modPSSR = (stat_counter / steps) % gl.config.PSSRatio;
  
  if (
      (gl.config.initialConditionT == 2)
      && (gl.config.evolveT == 1)
     )
    { heterogeneousQuench(gl, stat_counter, steps, originpoints, modSteps, modPSSR); } // heterogenous quench, hot blob T-profile 
  else
    { homogenousQuench(gl, stat_counter, steps, originpoints, modSteps, modPSSR); } // homogenous quench block
}

} // orch namespace ends here

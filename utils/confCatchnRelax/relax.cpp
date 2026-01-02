#define USE_MPI 
// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <string>
// #include <assert.h>

#include "plumbing/hila.h"

#include "glsol.hpp"
#include "matep.hpp"


void glsol::relax() {

  if (t < config.tdis && config.gamma.squarenorm() > 0 )
    {
      //hila::out0 << "config.gamma is " << config.gamma << "\n" << std::endl;

      /* this 2.0 infront gamma is not right */
      /* the despative term indeeds gives 2.0, however, this 2.0 later be absorbed into newly defined \gamma */
      pi[ALL] = pi[X] + (deltaPi[X] - 1.0 * config.gamma * pi[X])*config.dt; //Complex<real_t> C(a, b) = r + I *
      t += config.dt;
    }
  else
    {
      pi[ALL] = pi[X] + deltaPi[X]*config.dt;
      t += config.dt;
    }

} // relax() ends here


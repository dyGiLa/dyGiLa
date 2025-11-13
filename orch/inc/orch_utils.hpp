#ifndef ORCH_UTILS_HPP
#define ORCH_UTILS_HPP

//#define USE_PARIO 
#define USE_MPI 
// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
#include <string>
// #include <assert.h>

#include "plumbing/hila.h"
// #include "plumbing/globals.h" 

#include "glsol.hpp"
//#include "matep_namespace_utils.hpp" 

//#if defined USE_PARIO 
//#include "pario.hpp"
//#endif

namespace orch {
  void heterogeneousQuench(glsol &, unsigned int &, const unsigned int &, const CoordinateVector &, const unsigned int &, const unsigned int &);
  void homogenousQuench(glsol &, unsigned int &, const unsigned int &, const CoordinateVector &, const unsigned int &, const unsigned int &);  
} // namespace orch ends here

#endif

#ifndef DYGILA_UTILS_HPP
#define DYGILA_UTILS_HPP

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

namespace dyGiLa {
  namespace utils {
    void heterogeneousQuench(glsol &, unsigned int &, const unsigned int &, const CoordinateVector &, const unsigned int &, const unsigned int &);
    void homogenousQuench(glsol &, unsigned int &, const unsigned int &, const CoordinateVector &, const unsigned int &, const unsigned int &);
    void phaseMarking(glsol &, unsigned int &, const unsigned int &);
  } // namespace utils ends here
} // namespace dyGiLa ends here

#endif

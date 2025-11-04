#ifndef ORCH_HPP
#define ORCH_HPP

#define USE_PARIO 
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
#include "matep_namespace_utils.hpp" 

#if defined USE_PARIO 
#include "pario.hpp"
#endif

namespace orch {
  void writeHDF5_xmls(glsol &, parIO &);
  void pStreaming(glsol &, parIO &, unsigned int &, const unsigned int &);
  void gammaEvolve(glsol &, unsigned int &, const unsigned int &, const CoordinateVector &);
  void nextBlocks(glsol &, unsigned int &, const unsigned int &, const CoordinateVector &);
  std::tuple<glsol *const, const std::vector<std::string>, const CoordinateVector *const, const unsigned int> dyGiLaInit(int &, char **&);  
} // namespace orch ends here

#endif

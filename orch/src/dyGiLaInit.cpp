//#define USE_PARIO 
#define USE_MPI 
// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
#include <string>
// #include <assert.h>

#include "plumbing/hila.h"
#include "plumbing/globals.h" 

#include "glsol.hpp"
#include "matep_namespace_utils.hpp"
#include "orch.hpp"

// #if defined USE_PARIO 
// #include "pario.hpp"
// #endif

namespace orch {
  
std::tuple<glsol *const, const std::vector<std::string>, const CoordinateVector *const, const unsigned int> dyGiLaInit(int &argc, char **&argv) {
  glsol gl;

  std::vector<std::string> name_files = gl.allocate("sim_params.txt", argc, argv);

  // host & device memory initialization, gpuMemcpHostToDevice under hood
  matep::init_wrapper_mp();

  // initialize Temperature field
  gl.initializeT();

  // initilize static H-field
  gl.initializeH();

  // initialize OP field
  gl.initialize();

  // coordinate origin
  std::initializer_list<int> coordsList {0,0,0};
  const CoordinateVector originpoints(coordsList);

  // number of steps between reduction streaming
  const unsigned int steps = (gl.config.tEnd - gl.config.tStats)
                             /(gl.config.dt * gl.config.nOutputs);

  // return initialization list
  return std::make_tuple(&gl, name_files, &originpoints, steps);
  
} /* dyGiLaInit() end here */

} // orch namespace ends here


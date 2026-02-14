// #define USE_BOTHSIDE_GHOSTS
// #define USE_ADGRZ
// #define USE_MPI 
// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <string>
// #include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/fft.h"
//#include "plumbing/memalloc.h"

//#include "glsol.hpp"
//#include "matep.hpp"
#include "pario.hpp"

#include "ascent.hpp"
#include "conduit_blueprint.hpp"

void parIO::containerReserve_lVec() {
   
  l1Container.reserve(latticeVolumeWithGhost);
  l2Container.reserve(latticeVolumeWithGhost);
  l3Container.reserve(latticeVolumeWithGhost); 
  
} // containerReserve_lVec() end here


// #define USE_BOTHSIDE_GHOSTS
// #define USE_ADGRZ
// #define USE_MPI 
// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <string>
// #include <assert.h>

// #include "plumbing/hila.h"
// //#include "plumbing/fft.h"
// #include "plumbing/memalloc.h"

//#include "glsol.hpp"
//#include "matep.hpp"
#include "pario.hpp"

#include "ascent.hpp"
#include "conduit_blueprint.hpp"

void parIO::containerReserve_massCurrent() {
  
  jm1Container.reserve(latticeVolumeWithGhost);
  jm2Container.reserve(latticeVolumeWithGhost);
  jm3Container.reserve(latticeVolumeWithGhost);

  // phaseExpModulusContainer.reserve(latticeVolumeWithGhost);
  // phaseExpAngleContainer.reserve(latticeVolumeWithGhost);
  // phaseExp2ReContainer.reserve(latticeVolumeWithGhost);
  // phaseExp2ImContainer.reserve(latticeVolumeWithGhost);                  
    
} // containerReserve_massCurrent() end here


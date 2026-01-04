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
// #include "plumbing/fft.h"
// #include "plumbing/memalloc.h"

//#include "glsol.hpp"
//#include "matep.hpp"
#include "pario.hpp"

#include "ascent.hpp"
#include "conduit_blueprint.hpp"

void parIO::containerReserve_spinCurrent() {

  js11Container.reserve(latticeVolumeWithGhost);
  js21Container.reserve(latticeVolumeWithGhost);
  js31Container.reserve(latticeVolumeWithGhost);

  js12Container.reserve(latticeVolumeWithGhost);
  js22Container.reserve(latticeVolumeWithGhost);
  js32Container.reserve(latticeVolumeWithGhost);

  js13Container.reserve(latticeVolumeWithGhost);
  js23Container.reserve(latticeVolumeWithGhost);
  js33Container.reserve(latticeVolumeWithGhost);
      
} // containerReserve_spinCurrent() end here


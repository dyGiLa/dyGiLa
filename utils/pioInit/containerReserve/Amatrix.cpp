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
#include "plumbing/memalloc.h"

//#include "glsol.hpp"
//#include "matep.hpp"
#include "pario.hpp"

#include "ascent.hpp"
#include "conduit_blueprint.hpp"

void parIO::containerReserve_Amatrix() {
   
  u11Container.reserve(latticeVolumeWithGhost); v11Container.reserve(latticeVolumeWithGhost);
  u12Container.reserve(latticeVolumeWithGhost); v12Container.reserve(latticeVolumeWithGhost);
  u13Container.reserve(latticeVolumeWithGhost); v13Container.reserve(latticeVolumeWithGhost);
  u21Container.reserve(latticeVolumeWithGhost); v21Container.reserve(latticeVolumeWithGhost);
  u22Container.reserve(latticeVolumeWithGhost); v22Container.reserve(latticeVolumeWithGhost);
  u23Container.reserve(latticeVolumeWithGhost); v23Container.reserve(latticeVolumeWithGhost);
  u31Container.reserve(latticeVolumeWithGhost); v31Container.reserve(latticeVolumeWithGhost);
  u32Container.reserve(latticeVolumeWithGhost); v32Container.reserve(latticeVolumeWithGhost);
  u33Container.reserve(latticeVolumeWithGhost); v33Container.reserve(latticeVolumeWithGhost);    
         
} // containerReserve_Amatrix() end here


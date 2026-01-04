// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <string>
// #include <assert.h>

// #include "plumbing/hila.h"
// //#include "plumbing/fft.h"
// #include "plumbing/memalloc.h"

#include "glsol.hpp"
//#include "matep.hpp"
#include "pario.hpp"

#include "ascent.hpp"
#include "conduit_blueprint.hpp"

void parIO::containerReserve_gapAFETemPMarker(glsol &sol) {

    gapAContainer.reserve(latticeVolumeWithGhost);
    
    if (sol.config.pario_compute_feDensity == 1) { feDensityContainer.reserve(latticeVolumeWithGhost); }
    
    Temperature.reserve(latticeVolumeWithGhost);
    phaseMarker.reserve(latticeVolumeWithGhost);
    
} // containerReserve_gapAFETemPMarker() end here


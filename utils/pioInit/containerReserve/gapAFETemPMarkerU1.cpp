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

void parIO::containerReserve_gapAFETemPMarkerU1(glsol &sol) {

    gapAContainer.reserve(latticeVolumeWithGhost);
    
    if (sol.config.pario_compute_feDensity == 1) { feDensityContainer.reserve(latticeVolumeWithGhost); }    
    if (sol.config.pario_Temperature_pStream == 1) { Temperature.reserve(latticeVolumeWithGhost); }
    
    phaseMarker.reserve(latticeVolumeWithGhost);

    if (sol.config.pario_compute_U1Phase == 1) { U1_3phiContainer.reserve(latticeVolumeWithGhost); }
    
} // containerReserve_gapAFETemPMarkerU1() end here


#define USE_PARIO
#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
//#include <math.h>
#include <assert.h>

#include "plumbing/hila.h"
#include "plumbing/fft.h"

#include "glsol.hpp"
#include "matep.hpp"
#include "pario.hpp"

#include "ascent.hpp"
#include "conduit_blueprint.hpp"


void parIO::describeMesh_gapA_FEDensity() {

    // create an vertex associated field named gapAOrdered
    mesh["fields/gapAOrdered/association"] = "vertex";
    mesh["fields/gapAOrdered/topology"] = "topo";
    mesh["fields/gapAOrdered/values"].set_external(gapAOrdered.data(), latticeVolumeWithGhost);

    // create an vertex associated field named feDensityOrdered
    mesh["fields/feDensityOrdered/association"] = "vertex";
    mesh["fields/feDensityOrdered/topology"] = "topo";
    mesh["fields/feDensityOrdered/values"].set_external(feDensityOrdered.data(), latticeVolumeWithGhost);

} // describeMesh() end here


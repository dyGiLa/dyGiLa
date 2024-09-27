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
    mesh["fields/gapA/association"] = "vertex";
    mesh["fields/gapA/topology"] = "topo";
    mesh["fields/gapA/values"].set_external(gapAOrdered.data(), latticeVolumeWithGhost);

    // create an vertex associated field named feDensityOrdered
    mesh["fields/feDensity/association"] = "vertex";
    mesh["fields/feDensity/topology"] = "topo";
    mesh["fields/feDensity/values"].set_external(feDensityOrdered.data(), latticeVolumeWithGhost);

} // describeMesh() end here


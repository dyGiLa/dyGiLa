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


void parIO::describeMesh_phaseMarker() {

    // create an vertex associated field named phaseMarker
    mesh["fields/phaseMarker/association"] = "vertex";
    mesh["fields/phaseMarker/topology"] = "topo";
    mesh["fields/phaseMarker/values"].set_external(phaseMarker.data(), latticeVolumeWithGhost);

} // describeMesh() end here


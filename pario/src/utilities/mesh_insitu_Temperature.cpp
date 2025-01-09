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


void parIO::describeMesh_Temperature() {

    // create an vertex associated field named Temperature
    mesh["fields/Temperature/association"] = "vertex";
    mesh["fields/Temperature/topology"] = "topo";
    mesh["fields/Temperature/values"].set_external(Temperature.data(), latticeVolumeWithGhost);

} // describeMesh() end here


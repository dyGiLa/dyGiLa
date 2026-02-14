#define USE_MPI 
// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <string>
// #include <assert.h>

//#include "plumbing/hila.h"
//#include "plumbing/fft.h"

#include "glsol.hpp"
//#include "matep.hpp"
#include "pario.hpp"

#include "ascent.hpp"
#include "conduit_blueprint.hpp"


void parIO::describeMesh_lVec() {

    // create vertex associated field l1, l2, l3
    mesh["fields/l1/association"] = "vertex";
    mesh["fields/l1/topology"] = "topo";
    mesh["fields/l1/values"].set_external(l1Container.data(), latticeVolumeWithGhost);

    mesh["fields/l2/association"] = "vertex";
    mesh["fields/l2/topology"] = "topo";
    mesh["fields/l2/values"].set_external(l2Container.data(), latticeVolumeWithGhost);

    mesh["fields/l3/association"] = "vertex";
    mesh["fields/l3/topology"] = "topo";
    mesh["fields/l3/values"].set_external(l3Container.data(), latticeVolumeWithGhost);
        
} // describeMesh_lVector() end here


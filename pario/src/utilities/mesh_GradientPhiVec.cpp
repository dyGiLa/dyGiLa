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


void parIO::describeMesh_GradientPhiVec() {

    // create vertex associated field l1, l2, l3
    mesh["fields/GPhi1/association"] = "vertex";
    mesh["fields/GPhi1/topology"] = "topo";
    mesh["fields/GPhi1/values"].set_external(GPhi1Container.data(), latticeVolumeWithGhost);

    mesh["fields/GPhi2/association"] = "vertex";
    mesh["fields/GPhi2/topology"] = "topo";
    mesh["fields/GPhi2/values"].set_external(GPhi2Container.data(), latticeVolumeWithGhost);

    mesh["fields/GPhi3/association"] = "vertex";
    mesh["fields/GPhi3/topology"] = "topo";
    mesh["fields/GPhi3/values"].set_external(GPhi3Container.data(), latticeVolumeWithGhost);
        
} // describeMesh_GradientPhiVector() end here


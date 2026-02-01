#define USE_MPI 
// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <string>
#include <assert.h>

//#include "plumbing/hila.h"
//#include "plumbing/fft.h"

//#include "glsol.hpp"
//#include "matep.hpp"
#include "pario.hpp"

#include "ascent.hpp"
#include "conduit_blueprint.hpp"


void parIO::describeMesh_U13phi() {

    // create an vertex associated field named U1_3phi
    mesh["fields/U1_3phi/association"] = "vertex";
    mesh["fields/U1_3phi/topology"] = "topo";
    mesh["fields/U1_3phi/values"].set_external(U1_3phiContainer.data(), latticeVolumeWithGhost);

} // describeMesh_U13phi() end here


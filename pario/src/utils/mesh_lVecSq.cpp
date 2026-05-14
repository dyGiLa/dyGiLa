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


void parIO::describeMesh_lVecSq() {

    mesh["fields/l_Sq/association"] = "vertex";
    mesh["fields/l_Sq/topology"] = "topo";
    mesh["fields/l_Sq/values"].set_external(lsqContainer.data(), latticeVolumeWithGhost);
        
} // describeMesh_lVector() end here


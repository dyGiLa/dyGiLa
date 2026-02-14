#define USE_MPI 
// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <string>
// #include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/fft.h"

//#include "glsol.hpp"
//#include "matep.hpp"
#include "pario.hpp"

#include "ascent.hpp"
#include "conduit_blueprint.hpp"


void parIO::describeMesh_massCurrent() {

      mesh["fields/jm1Container/association"] = "vertex";
      mesh["fields/jm1Container/topology"] = "topo";
      mesh["fields/jm1Container/values"].set_external(jm1Container.data(), latticeVolumeWithGhost);

      mesh["fields/jm2Container/association"] = "vertex";
      mesh["fields/jm2Container/topology"] = "topo";
      mesh["fields/jm2Container/values"].set_external(jm2Container.data(), latticeVolumeWithGhost);

      mesh["fields/jm3Container/association"] = "vertex";
      mesh["fields/jm3Container/topology"] = "topo";
      mesh["fields/jm3Container/values"].set_external(jm3Container.data(), latticeVolumeWithGhost);

} // describeMesh() end here


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


void parIO::describeMesh_massCurrent() {

      mesh["fields/jm1Ordered/association"] = "vertex";
      mesh["fields/jm1Ordered/topology"] = "topo";
      mesh["fields/jm1Ordered/values"].set_external(jm1Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/jm2Ordered/association"] = "vertex";
      mesh["fields/jm2Ordered/topology"] = "topo";
      mesh["fields/jm2Ordered/values"].set_external(jm2Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/jm3Ordered/association"] = "vertex";
      mesh["fields/jm3Ordered/topology"] = "topo";
      mesh["fields/jm3Ordered/values"].set_external(jm3Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/phaseExpModulusO/association"] = "vertex";
      mesh["fields/phaseExpModulusO/topology"] = "topo";
      mesh["fields/phaseExpModulusO/values"].set_external(phaseExpModulusO.data(), latticeVolumeWithGhost);

      mesh["fields/phaseExpAngleO/association"] = "vertex";
      mesh["fields/phaseExpAngleO/topology"] = "topo";
      mesh["fields/phaseExpAngleO/values"].set_external(phaseExpAngleO.data(), latticeVolumeWithGhost);

      mesh["fields/phaseExp2ReO/association"] = "vertex";
      mesh["fields/phaseExp2ReO/topology"] = "topo";
      mesh["fields/phaseExp2ReO/values"].set_external(phaseExp2ReO.data(), latticeVolumeWithGhost);

      mesh["fields/phaseExp2ImO/association"] = "vertex";
      mesh["fields/phaseExp2ImO/topology"] = "topo";
      mesh["fields/phaseExp2ImO/values"].set_external(phaseExp2ImO.data(), latticeVolumeWithGhost);      

} // describeMesh() end here


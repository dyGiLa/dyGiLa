#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>

#include "plumbing/hila.h"
#include "plumbing/fft.h"

#include "glsol.hpp"
#include "matep.hpp"
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

      mesh["fields/phaseExpModulusContainer/association"] = "vertex";
      mesh["fields/phaseExpModulusContainer/topology"] = "topo";
      mesh["fields/phaseExpModulusContainer/values"].set_external(phaseExpModulusContainer.data(), latticeVolumeWithGhost);

      mesh["fields/phaseExpAngleContainer/association"] = "vertex";
      mesh["fields/phaseExpAngleContainer/topology"] = "topo";
      mesh["fields/phaseExpAngleContainer/values"].set_external(phaseExpAngleContainer.data(), latticeVolumeWithGhost);

      mesh["fields/phaseExp2ReContainer/association"] = "vertex";
      mesh["fields/phaseExp2ReContainer/topology"] = "topo";
      mesh["fields/phaseExp2ReContainer/values"].set_external(phaseExp2ReContainer.data(), latticeVolumeWithGhost);

      mesh["fields/phaseExp2ImContainer/association"] = "vertex";
      mesh["fields/phaseExp2ImContainer/topology"] = "topo";
      mesh["fields/phaseExp2ImContainer/values"].set_external(phaseExp2ImContainer.data(), latticeVolumeWithGhost);      

} // describeMesh() end here


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


void parIO::describeMesh_spinCurrent() {

      mesh["fields/js11Container/association"] = "vertex";
      mesh["fields/js11Container/topology"] = "topo";
      mesh["fields/js11Container/values"].set_external(js11Container.data(), latticeVolumeWithGhost);

      mesh["fields/js21Container/association"] = "vertex";
      mesh["fields/js21Container/topology"] = "topo";
      mesh["fields/js21Container/values"].set_external(js21Container.data(), latticeVolumeWithGhost);

      mesh["fields/js31Container/association"] = "vertex";
      mesh["fields/js31Container/topology"] = "topo";
      mesh["fields/js31Container/values"].set_external(js31Container.data(), latticeVolumeWithGhost);

      mesh["fields/js12Container/association"] = "vertex";
      mesh["fields/js12Container/topology"] = "topo";
      mesh["fields/js12Container/values"].set_external(js12Container.data(), latticeVolumeWithGhost);

      mesh["fields/js22Container/association"] = "vertex";
      mesh["fields/js22Container/topology"] = "topo";
      mesh["fields/js22Container/values"].set_external(js22Container.data(), latticeVolumeWithGhost);

      mesh["fields/js32Container/association"] = "vertex";
      mesh["fields/js32Container/topology"] = "topo";
      mesh["fields/js32Container/values"].set_external(js32Container.data(), latticeVolumeWithGhost);

      mesh["fields/js13Container/association"] = "vertex";
      mesh["fields/js13Container/topology"] = "topo";
      mesh["fields/js13Container/values"].set_external(js13Container.data(), latticeVolumeWithGhost);

      mesh["fields/js23Container/association"] = "vertex";
      mesh["fields/js23Container/topology"] = "topo";
      mesh["fields/js23Container/values"].set_external(js23Container.data(), latticeVolumeWithGhost);

      mesh["fields/js33Container/association"] = "vertex";
      mesh["fields/js33Container/topology"] = "topo";
      mesh["fields/js33Container/values"].set_external(js33Container.data(), latticeVolumeWithGhost);
   
} // describeMesh() end here


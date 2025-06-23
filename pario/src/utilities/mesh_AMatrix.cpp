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


void parIO::describeMesh_AMatrix() {

    // /*----------------------------------------------------------------------*/
    // /*---create vertices associated field named of uxxOrdered vxxOrdered ---*/
    // /*----------------------------------------------------------------------*/
      mesh["fields/u11Container/association"] = "vertex";
      mesh["fields/u11Container/topology"] = "topo";
      mesh["fields/u11Container/values"].set_external(u11Container.data(), latticeVolumeWithGhost);

      mesh["fields/u12Container/association"] = "vertex";
      mesh["fields/u12Container/topology"] = "topo";
      mesh["fields/u12Container/values"].set_external(u12Container.data(), latticeVolumeWithGhost);

      mesh["fields/u13Container/association"] = "vertex";
      mesh["fields/u13Container/topology"] = "topo";
      mesh["fields/u13Container/values"].set_external(u13Container.data(), latticeVolumeWithGhost);

      mesh["fields/u21Container/association"] = "vertex";
      mesh["fields/u21Container/topology"] = "topo";
      mesh["fields/u21Container/values"].set_external(u21Container.data(), latticeVolumeWithGhost);

      mesh["fields/u22Container/association"] = "vertex";
      mesh["fields/u22Container/topology"] = "topo";
      mesh["fields/u22Container/values"].set_external(u22Container.data(), latticeVolumeWithGhost);

      mesh["fields/u23Container/association"] = "vertex";
      mesh["fields/u23Container/topology"] = "topo";
      mesh["fields/u23Container/values"].set_external(u23Container.data(), latticeVolumeWithGhost);

      mesh["fields/u31Container/association"] = "vertex";
      mesh["fields/u31Container/topology"] = "topo";
      mesh["fields/u31Container/values"].set_external(u31Container.data(), latticeVolumeWithGhost);

      mesh["fields/u32Container/association"] = "vertex";
      mesh["fields/u32Container/topology"] = "topo";
      mesh["fields/u32Container/values"].set_external(u32Container.data(), latticeVolumeWithGhost);

      mesh["fields/u33Container/association"] = "vertex";
      mesh["fields/u33Container/topology"] = "topo";
      mesh["fields/u33Container/values"].set_external(u33Container.data(), latticeVolumeWithGhost);

      mesh["fields/v11Container/association"] = "vertex";
      mesh["fields/v11Container/topology"] = "topo";
      mesh["fields/v11Container/values"].set_external(v11Container.data(), latticeVolumeWithGhost);

      mesh["fields/v12Container/association"] = "vertex";
      mesh["fields/v12Container/topology"] = "topo";
      mesh["fields/v12Container/values"].set_external(v12Container.data(), latticeVolumeWithGhost);

      mesh["fields/v13Container/association"] = "vertex";
      mesh["fields/v13Container/topology"] = "topo";
      mesh["fields/v13Container/values"].set_external(v13Container.data(), latticeVolumeWithGhost);

      mesh["fields/v21Container/association"] = "vertex";
      mesh["fields/v21Container/topology"] = "topo";
      mesh["fields/v21Container/values"].set_external(v21Container.data(), latticeVolumeWithGhost);

      mesh["fields/v22Container/association"] = "vertex";
      mesh["fields/v22Container/topology"] = "topo";
      mesh["fields/v22Container/values"].set_external(v22Container.data(), latticeVolumeWithGhost);

      mesh["fields/v23Container/association"] = "vertex";
      mesh["fields/v23Container/topology"] = "topo";
      mesh["fields/v23Container/values"].set_external(v23Container.data(), latticeVolumeWithGhost);

      mesh["fields/v31Container/association"] = "vertex";
      mesh["fields/v31Container/topology"] = "topo";
      mesh["fields/v31Container/values"].set_external(v31Container.data(), latticeVolumeWithGhost);

      mesh["fields/v32Container/association"] = "vertex";
      mesh["fields/v32Container/topology"] = "topo";
      mesh["fields/v32Container/values"].set_external(v32Container.data(), latticeVolumeWithGhost);

      mesh["fields/v33Container/association"] = "vertex";
      mesh["fields/v33Container/topology"] = "topo";
      mesh["fields/v33Container/values"].set_external(v33Container.data(), latticeVolumeWithGhost);              

} // describeMesh() end here


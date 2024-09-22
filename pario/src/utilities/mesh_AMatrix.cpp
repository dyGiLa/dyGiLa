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


void parIO::describeMesh_AMatrix() {

    // /*----------------------------------------------------------------------*/
    // /*---create vertices associated field named of uxxOrdered vxxOrdered ---*/
    // /*----------------------------------------------------------------------*/
      mesh["fields/u11Ordered/association"] = "vertex";
      mesh["fields/u11Ordered/topology"] = "topo";
      mesh["fields/u11Ordered/values"].set_external(u11Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/u12Ordered/association"] = "vertex";
      mesh["fields/u12Ordered/topology"] = "topo";
      mesh["fields/u12Ordered/values"].set_external(u12Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/u13Ordered/association"] = "vertex";
      mesh["fields/u13Ordered/topology"] = "topo";
      mesh["fields/u13Ordered/values"].set_external(u13Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/u21Ordered/association"] = "vertex";
      mesh["fields/u21Ordered/topology"] = "topo";
      mesh["fields/u21Ordered/values"].set_external(u21Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/u22Ordered/association"] = "vertex";
      mesh["fields/u22Ordered/topology"] = "topo";
      mesh["fields/u22Ordered/values"].set_external(u22Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/u23Ordered/association"] = "vertex";
      mesh["fields/u23Ordered/topology"] = "topo";
      mesh["fields/u23Ordered/values"].set_external(u23Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/u31Ordered/association"] = "vertex";
      mesh["fields/u31Ordered/topology"] = "topo";
      mesh["fields/u31Ordered/values"].set_external(u31Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/u32Ordered/association"] = "vertex";
      mesh["fields/u32Ordered/topology"] = "topo";
      mesh["fields/u32Ordered/values"].set_external(u32Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/u33Ordered/association"] = "vertex";
      mesh["fields/u33Ordered/topology"] = "topo";
      mesh["fields/u33Ordered/values"].set_external(u33Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/v11Ordered/association"] = "vertex";
      mesh["fields/v11Ordered/topology"] = "topo";
      mesh["fields/v11Ordered/values"].set_external(v11Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/v12Ordered/association"] = "vertex";
      mesh["fields/v12Ordered/topology"] = "topo";
      mesh["fields/v12Ordered/values"].set_external(v12Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/v13Ordered/association"] = "vertex";
      mesh["fields/v13Ordered/topology"] = "topo";
      mesh["fields/v13Ordered/values"].set_external(v13Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/v21Ordered/association"] = "vertex";
      mesh["fields/v21Ordered/topology"] = "topo";
      mesh["fields/v21Ordered/values"].set_external(v21Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/v22Ordered/association"] = "vertex";
      mesh["fields/v22Ordered/topology"] = "topo";
      mesh["fields/v22Ordered/values"].set_external(v22Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/v23Ordered/association"] = "vertex";
      mesh["fields/v23Ordered/topology"] = "topo";
      mesh["fields/v23Ordered/values"].set_external(v23Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/v31Ordered/association"] = "vertex";
      mesh["fields/v31Ordered/topology"] = "topo";
      mesh["fields/v31Ordered/values"].set_external(v31Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/v32Ordered/association"] = "vertex";
      mesh["fields/v32Ordered/topology"] = "topo";
      mesh["fields/v32Ordered/values"].set_external(v32Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/v33Ordered/association"] = "vertex";
      mesh["fields/v33Ordered/topology"] = "topo";
      mesh["fields/v33Ordered/values"].set_external(v33Ordered.data(), latticeVolumeWithGhost);              

} // describeMesh() end here


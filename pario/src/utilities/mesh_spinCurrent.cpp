#define USE_PARIO
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


void parIO::describeMesh_spinCurrent() {

      mesh["fields/js11O/association"] = "vertex";
      mesh["fields/js11O/topology"] = "topo";
      mesh["fields/js11O/values"].set_external(js11O.data(), latticeVolumeWithGhost);

      mesh["fields/js21O/association"] = "vertex";
      mesh["fields/js21O/topology"] = "topo";
      mesh["fields/js21O/values"].set_external(js21O.data(), latticeVolumeWithGhost);

      mesh["fields/js31O/association"] = "vertex";
      mesh["fields/js31O/topology"] = "topo";
      mesh["fields/js31O/values"].set_external(js31O.data(), latticeVolumeWithGhost);

      mesh["fields/js12O/association"] = "vertex";
      mesh["fields/js12O/topology"] = "topo";
      mesh["fields/js12O/values"].set_external(js12O.data(), latticeVolumeWithGhost);

      mesh["fields/js22O/association"] = "vertex";
      mesh["fields/js22O/topology"] = "topo";
      mesh["fields/js22O/values"].set_external(js22O.data(), latticeVolumeWithGhost);

      mesh["fields/js32O/association"] = "vertex";
      mesh["fields/js32O/topology"] = "topo";
      mesh["fields/js32O/values"].set_external(js32O.data(), latticeVolumeWithGhost);

      mesh["fields/js13O/association"] = "vertex";
      mesh["fields/js13O/topology"] = "topo";
      mesh["fields/js13O/values"].set_external(js13O.data(), latticeVolumeWithGhost);

      mesh["fields/js23O/association"] = "vertex";
      mesh["fields/js23O/topology"] = "topo";
      mesh["fields/js23O/values"].set_external(js23O.data(), latticeVolumeWithGhost);

      mesh["fields/js33O/association"] = "vertex";
      mesh["fields/js33O/topology"] = "topo";
      mesh["fields/js33O/values"].set_external(js33O.data(), latticeVolumeWithGhost);
   
} // describeMesh() end here


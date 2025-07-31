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


void parIO::defineActions_spinCurrent(glsol &sol) {

      conduit::Node &add_act7 = actions.append();
      add_act7["action"] = "add_extracts";

      conduit::Node &extracts = add_act7["extracts"];
      extracts["e3/type"] = "relay";
      extracts["e3/params/path"] = "pio_Current/dyGiLa-sim-spinCurrent";
      extracts["e3/params/protocol"] = "blueprint/mesh/hdf5";

      // extracts["e3/params/fields"].append().set("gapAContainer");
      // extracts["e3/params/fields"].append().set("feDensityContainer");
      extracts["e3/params/fields"].append().set("gapA");
      if (sol.config.pario_compute_feDensity == 1) { extracts["e3/params/fields"].append().set("feDensityr"); }
      
      extracts["e3/params/fields"].append().set("js11Container");
      extracts["e3/params/fields"].append().set("js21Container");
      extracts["e3/params/fields"].append().set("js31Container");
      extracts["e3/params/fields"].append().set("js12Container");
      extracts["e3/params/fields"].append().set("js22Container");
      extracts["e3/params/fields"].append().set("js32Container");
      extracts["e3/params/fields"].append().set("js13Container");
      extracts["e3/params/fields"].append().set("js23Container");
      extracts["e3/params/fields"].append().set("js33Container");      
    
} // defineActions() call end here


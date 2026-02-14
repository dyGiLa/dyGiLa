#define USE_MPI 
// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <string>
// #include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/fft.h"

#include "glsol.hpp"
//#include "matep.hpp"
#include "pario.hpp"

#include "ascent.hpp"
#include "conduit_blueprint.hpp"


void parIO::defineActions_massCurrent(glsol &sol) {

      conduit::Node &add_act6 = actions.append();
      add_act6["action"] = "add_extracts";

      conduit::Node &extracts = add_act6["extracts"];
      extracts["e2/type"] = "relay";
      extracts["e2/params/path"] = "pio_Current/dyGiLa-sim-massCurrent";
      extracts["e2/params/protocol"] = "blueprint/mesh/hdf5";

      // extracts["e2/params/fields"].append().set("gapAContainer");
      // extracts["e2/params/fields"].append().set("feDensityContainer");
      extracts["e2/params/fields"].append().set("gapA");
      if (sol.config.pario_compute_feDensity == 1) { extracts["e2/params/fields"].append().set("feDensity"); }
      
      extracts["e2/params/fields"].append().set("jm1Container");
      extracts["e2/params/fields"].append().set("jm2Container");
      extracts["e2/params/fields"].append().set("jm3Container");
    
} // defineActions() call end here


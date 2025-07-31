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


void parIO::defineActions_phaseMarker(glsol &sol) {

    // /* >>>>>>>>>>>>> extract hdf5 <<<<<<<<<<<<<< */

     conduit::Node &add_act8 = actions.append();
     add_act8["action"] = "add_extracts";

     conduit::Node &extracts = add_act8["extracts"];
     extracts["e4/type"] = "relay";
     extracts["e4/params/path"] = "pio/dyGiLa-sim-pMarker";
     extracts["e4/params/protocol"] = "blueprint/mesh/hdf5";

     extracts["e4/params/fields"].append().set("gapA");
     if (sol.config.pario_compute_feDensity == 1) { extracts["e4/params/fields"].append().set("feDensity"); }
     
     extracts["e4/params/fields"].append().set("phaseMarker");
    
} // defineActions() call end here


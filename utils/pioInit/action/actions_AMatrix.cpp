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


void parIO::defineActions_AMatrix(glsol &sol) {

    // /* >>>>>>>>>>>>> extract hdf5 <<<<<<<<<<<<<< */

     conduit::Node &add_act5 = actions.append();
     add_act5["action"] = "add_extracts";

     conduit::Node &extracts = add_act5["extracts"];
     extracts["e1/type"] = "relay";
     extracts["e1/params/path"] = "pio/dyGiLa-sim-Amatrix_t-%09d";
     extracts["e1/params/protocol"] = "blueprint/mesh/hdf5";
     // extracts["e1/params/protocol"] = "hdf5";     

     extracts["e1/params/fields"].append().set("gapA");
     if (sol.config.pario_compute_feDensity == 1) { extracts["e1/params/fields"].append().set("feDensity"); }
    
     extracts["e1/params/fields"].append().set("u11Container");
     extracts["e1/params/fields"].append().set("u12Container");
     extracts["e1/params/fields"].append().set("u13Container");
     extracts["e1/params/fields"].append().set("u21Container");
     extracts["e1/params/fields"].append().set("u22Container");
     extracts["e1/params/fields"].append().set("u23Container");
     extracts["e1/params/fields"].append().set("u31Container");
     extracts["e1/params/fields"].append().set("u32Container");
     extracts["e1/params/fields"].append().set("u33Container");

     extracts["e1/params/fields"].append().set("v11Container");
     extracts["e1/params/fields"].append().set("v12Container");
     extracts["e1/params/fields"].append().set("v13Container");
     extracts["e1/params/fields"].append().set("v21Container");
     extracts["e1/params/fields"].append().set("v22Container");
     extracts["e1/params/fields"].append().set("v23Container");
     extracts["e1/params/fields"].append().set("v31Container");
     extracts["e1/params/fields"].append().set("v32Container");
     extracts["e1/params/fields"].append().set("v33Container");
     
} // defineActions() call end here


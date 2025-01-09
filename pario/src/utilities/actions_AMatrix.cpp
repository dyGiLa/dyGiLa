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


void parIO::defineActions_AMatrix() {

    // /* >>>>>>>>>>>>> extract hdf5 <<<<<<<<<<<<<< */

     conduit::Node &add_act3 = actions.append();
     add_act3["action"] = "add_extracts";

     conduit::Node &extracts = add_act3["extracts"];
     extracts["e1/type"] = "relay";
     extracts["e1/params/path"] = "dyGiLa-sim-data";
     extracts["e1/params/protocol"] = "blueprint/mesh/hdf5";

     // extracts["e1/params/fields"].append().set("gapA");
     // extracts["e1/params/fields"].append().set("feDensity");
    
     extracts["e1/params/fields"].append().set("u11Ordered");
     extracts["e1/params/fields"].append().set("u12Ordered");
     extracts["e1/params/fields"].append().set("u13Ordered");
     extracts["e1/params/fields"].append().set("u21Ordered");
     extracts["e1/params/fields"].append().set("u22Ordered");
     extracts["e1/params/fields"].append().set("u23Ordered");
     extracts["e1/params/fields"].append().set("u31Ordered");
     extracts["e1/params/fields"].append().set("u32Ordered");
     extracts["e1/params/fields"].append().set("u33Ordered");

     extracts["e1/params/fields"].append().set("v11Ordered");
     extracts["e1/params/fields"].append().set("v12Ordered");
     extracts["e1/params/fields"].append().set("v13Ordered");
     extracts["e1/params/fields"].append().set("v21Ordered");
     extracts["e1/params/fields"].append().set("v22Ordered");
     extracts["e1/params/fields"].append().set("v23Ordered");
     extracts["e1/params/fields"].append().set("v31Ordered");
     extracts["e1/params/fields"].append().set("v32Ordered");
     extracts["e1/params/fields"].append().set("v33Ordered");
    
} // defineActions() call end here


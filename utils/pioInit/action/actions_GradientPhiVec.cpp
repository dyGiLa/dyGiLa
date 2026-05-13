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


void parIO::defineActions_GradientPhiVec() {

      /* this is fuction only for relay GPhi extact, not for insitu */
  
      conduit::Node &add_act10 = actions.append();
      add_act10["action"] = "add_extracts";

      conduit::Node &extracts = add_act10["extracts"];
      extracts["e6/type"] = "relay";
      extracts["e6/params/path"] = "pio_Vec/dyGiLa-sim-GPhiVector_t-%09d";
      extracts["e6/params/protocol"] = "blueprint/mesh/hdf5";

      extracts["e6/params/fields"].append().set("GPhi1");
      extracts["e6/params/fields"].append().set("GPhi2");
      extracts["e6/params/fields"].append().set("GPhi3");
    
} // defineActions() call end here


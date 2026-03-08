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


void parIO::defineActions_lVec() {

      /* this is fuction only for relay l_i extact, not for insitu */
  
      conduit::Node &add_act9 = actions.append();
      add_act9["action"] = "add_extracts";

      conduit::Node &extracts = add_act9["extracts"];
      extracts["e5/type"] = "relay";
      extracts["e5/params/path"] = "pio_Vec/dyGiLa-sim-lVector_t-%09d";
      extracts["e5/params/protocol"] = "blueprint/mesh/hdf5";

      //extracts["e5/params/fields"].append().set("l_Sq");
      extracts["e5/params/fields"].append().set("l1");
      extracts["e5/params/fields"].append().set("l2");
      extracts["e5/params/fields"].append().set("l3");
    
} // defineActions() call end here


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


void parIO::defineActions_spinCurrent() {

      conduit::Node &add_act7 = actions.append();
      add_act7["action"] = "add_extracts";

      conduit::Node &extracts = add_act7["extracts"];
      extracts["e1/type"] = "relay";
      extracts["e1/params/path"] = "sim-data";
      extracts["e1/params/protocol"] = "blueprint/mesh/hdf5";

      extracts["e1/params/fields"].append().set("gapAContainer");
      extracts["e1/params/fields"].append().set("feDensityContainer");
      extracts["e1/params/fields"].append().set("js11Container");
      extracts["e1/params/fields"].append().set("js21Container");
      extracts["e1/params/fields"].append().set("js31Container");
      extracts["e1/params/fields"].append().set("js12Container");
      extracts["e1/params/fields"].append().set("js22Container");
      extracts["e1/params/fields"].append().set("js32Container");
      extracts["e1/params/fields"].append().set("js13Container");
      extracts["e1/params/fields"].append().set("js23Container");
      extracts["e1/params/fields"].append().set("js33Container");      
    
} // defineActions() call end here


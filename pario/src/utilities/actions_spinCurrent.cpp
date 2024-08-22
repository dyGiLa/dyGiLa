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


void parIO::defineActions_spinCurrent() {

      conduit::Node &add_act7 = actions.append();
      add_act7["action"] = "add_extracts";

      conduit::Node &extracts = add_act7["extracts"];
      extracts["e1/type"] = "relay";
      extracts["e1/params/path"] = "sim-data";
      extracts["e1/params/protocol"] = "blueprint/mesh/hdf5";

      extracts["e1/params/fields"].append().set("gapAOrdered");
      extracts["e1/params/fields"].append().set("feDensityOrdered");
      extracts["e1/params/fields"].append().set("js11O");
      extracts["e1/params/fields"].append().set("js21O");
      extracts["e1/params/fields"].append().set("js31O");
      extracts["e1/params/fields"].append().set("js12O");
      extracts["e1/params/fields"].append().set("js22O");
      extracts["e1/params/fields"].append().set("js32O");
      extracts["e1/params/fields"].append().set("js13O");
      extracts["e1/params/fields"].append().set("js23O");
      extracts["e1/params/fields"].append().set("js33O");      
    
} // defineActions() call end here


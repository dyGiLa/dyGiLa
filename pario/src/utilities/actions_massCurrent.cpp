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


void parIO::defineActions_massCurrent() {

      conduit::Node &add_act6 = actions.append();
      add_act6["action"] = "add_extracts";

      conduit::Node &extracts = add_act6["extracts"];
      extracts["e1/type"] = "relay";
      extracts["e1/params/path"] = "sim-data";
      extracts["e1/params/protocol"] = "blueprint/mesh/hdf5";

      extracts["e1/params/fields"].append().set("gapAOrdered");
      extracts["e1/params/fields"].append().set("feDensityOrdered");
      extracts["e1/params/fields"].append().set("jm1Ordered");
      extracts["e1/params/fields"].append().set("jm2Ordered");
      extracts["e1/params/fields"].append().set("jm3Ordered");
      extracts["e1/params/fields"].append().set("phaseExpModulusO");
      extracts["e1/params/fields"].append().set("phaseExpAngleO");
      extracts["e1/params/fields"].append().set("phaseExp2ReO");
      extracts["e1/params/fields"].append().set("phaseExp2ImO");                       
    
} // defineActions() call end here


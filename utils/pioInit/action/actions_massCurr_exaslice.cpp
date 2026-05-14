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


void parIO::defineActions_massCurrent_exaslice(glsol &sol) {

      /* this is fuction only for relay extact, not for insitu */

      conduit::Node &add_act8 = actions.append();
      add_act8["action"] = "add_pipelines";
      conduit::Node &pipelines6 = add_act8["pipelines"];
      pipelines6["pl18/f1/type"] = "exaslice";
      conduit::Node &slice_params = pipelines6["pl18/f1/params"];
      slice_params["point/x"] = sol.config.massCurr_slice1_point_x;
      slice_params["point/y"] = sol.config.massCurr_slice1_point_y;
      slice_params["point/z"] = sol.config.massCurr_slice1_point_z;
      slice_params["normal/x"] = sol.config.massCurr_slice1_norm_x;
      slice_params["normal/y"] = sol.config.massCurr_slice1_norm_y;
      slice_params["normal/z"] = sol.config.massCurr_slice1_norm_z;
  
      conduit::Node &add_act14 = actions.append();
      add_act14["action"] = "add_extracts";
      conduit::Node &extracts = add_act14["extracts"];
      extracts["e10/type"] = "relay";
      extracts["e10/pipeline"]  = "pl18";
      extracts["e10/params/path"] = "pio/dyGiLa-sim-massCurr-exaslice_t-%09d";
      extracts["e10/params/protocol"] = "blueprint/mesh/hdf5";
      
      extracts["e10/params/fields"].append().set("jm1Container");
      extracts["e10/params/fields"].append().set("jm2Container");
      extracts["e10/params/fields"].append().set("jm3Container");
      
} // defineActions() call end here


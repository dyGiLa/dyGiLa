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


void parIO::defineActions_spinCurrent_exaslice(glsol &sol) {

      /* this is fuction only for relay extact, not for insitu */

      conduit::Node &add_act7 = actions.append();
      add_act7["action"] = "add_pipelines";
      conduit::Node &pipelines5 = add_act7["pipelines"];
      pipelines5["pl17/f1/type"] = "exaslice";
      conduit::Node &slice_params = pipelines5["pl17/f1/params"];
      slice_params["point/x"] = sol.config.spinCurr_slice1_point_x;
      slice_params["point/y"] = sol.config.spinCurr_slice1_point_y;
      slice_params["point/z"] = sol.config.spinCurr_slice1_point_z;
      slice_params["normal/x"] = sol.config.spinCurr_slice1_norm_x;
      slice_params["normal/y"] = sol.config.spinCurr_slice1_norm_y;
      slice_params["normal/z"] = sol.config.spinCurr_slice1_norm_z;
  
      conduit::Node &add_act13 = actions.append();
      add_act13["action"] = "add_extracts";
      conduit::Node &extracts = add_act13["extracts"];
      extracts["e9/type"] = "relay";
      extracts["e9/pipeline"]  = "pl17";
      extracts["e9/params/path"] = "pio/dyGiLa-sim-spinCurr-exaslice_t-%09d";
      extracts["e9/params/protocol"] = "blueprint/mesh/hdf5";
      
      extracts["e9/params/fields"].append().set("js11Container");
      extracts["e9/params/fields"].append().set("js21Container");
      extracts["e9/params/fields"].append().set("js31Container");
      extracts["e9/params/fields"].append().set("js12Container");
      extracts["e9/params/fields"].append().set("js22Container");
      extracts["e9/params/fields"].append().set("js32Container");
      extracts["e9/params/fields"].append().set("js13Container");
      extracts["e9/params/fields"].append().set("js23Container");
      extracts["e9/params/fields"].append().set("js33Container");      
      
} // defineActions() call end here


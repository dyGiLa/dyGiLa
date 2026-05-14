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


void parIO::defineActions_AMatrix_exaslice(glsol &sol) {

      /* this is fuction only for relay AM extact, not for insitu */

      conduit::Node &add_act6 = actions.append();
      add_act6["action"] = "add_pipelines";
      conduit::Node &pipelines4 = add_act6["pipelines"];
      pipelines4["pl16/f1/type"] = "exaslice";
      conduit::Node &slice_params = pipelines4["pl16/f1/params"];
      slice_params["point/x"] = sol.config.AM_slice1_point_x;
      slice_params["point/y"] = sol.config.AM_slice1_point_y;
      slice_params["point/z"] = sol.config.AM_slice1_point_z;
      slice_params["normal/x"] = sol.config.AM_slice1_norm_x;
      slice_params["normal/y"] = sol.config.AM_slice1_norm_y;
      slice_params["normal/z"] = sol.config.AM_slice1_norm_z;
  
      conduit::Node &add_act12 = actions.append();
      add_act12["action"] = "add_extracts";
      conduit::Node &extracts = add_act12["extracts"];
      extracts["e8/type"] = "relay";
      extracts["e8/pipeline"]  = "pl16";
      extracts["e8/params/path"] = "pio/dyGiLa-sim-AMatrix-exaslice_t-%09d";
      extracts["e8/params/protocol"] = "blueprint/mesh/hdf5";
      
      extracts["e8/params/fields"].append().set("u11Container");
      extracts["e8/params/fields"].append().set("u12Container");
      extracts["e8/params/fields"].append().set("u13Container");
      extracts["e8/params/fields"].append().set("u21Container");
      extracts["e8/params/fields"].append().set("u22Container");
      extracts["e8/params/fields"].append().set("u23Container");
      extracts["e8/params/fields"].append().set("u31Container");
      extracts["e8/params/fields"].append().set("u32Container");
      extracts["e8/params/fields"].append().set("u33Container");

      extracts["e8/params/fields"].append().set("v11Container");
      extracts["e8/params/fields"].append().set("v12Container");
      extracts["e8/params/fields"].append().set("v13Container");
      extracts["e8/params/fields"].append().set("v21Container");
      extracts["e8/params/fields"].append().set("v22Container");
      extracts["e8/params/fields"].append().set("v23Container");
      extracts["e8/params/fields"].append().set("v31Container");
      extracts["e8/params/fields"].append().set("v32Container");
      extracts["e8/params/fields"].append().set("v33Container");
          
} // defineActions() call end here


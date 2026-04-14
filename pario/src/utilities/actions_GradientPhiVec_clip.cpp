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


void parIO::defineActions_GradientPhiVec_clip(glsol &sol) {

      /* this is fuction only for relay GPhi extact, not for insitu */

      conduit::Node &add_act6 = actions.append();
      add_act6["action"] = "add_pipelines";
      conduit::Node &pipelines4 = add_act6["pipelines"];
      pipelines4["pl16/f1/type"] = "clip";
      conduit::Node &clip_params = pipelines4["pl16/f1/params"];
      clip_params["topology"] = "topo";
      clip_params["plane/point/x"] = sol.config.GPhi_clip1_point_x;
      clip_params["plane/point/y"] = sol.config.GPhi_clip1_point_y;
      clip_params["plane/point/z"] = sol.config.GPhi_clip1_point_z;
      clip_params["plane/normal/x"] = sol.config.GPhi_clip1_norm_x;
      clip_params["plane/normal/y"] = sol.config.GPhi_clip1_norm_y;
      clip_params["plane/normal/z"] = sol.config.GPhi_clip1_norm_z;
  
      conduit::Node &add_act12 = actions.append();
      add_act12["action"] = "add_extracts";
      conduit::Node &extracts = add_act12["extracts"];
      extracts["e8/type"] = "relay";
      extracts["e8/pipeline"]  = "pl16";
      extracts["e8/params/path"] = "pio_Vec/dyGiLa-sim-GPhiVector-clip_t-%09d";
      extracts["e8/params/protocol"] = "blueprint/mesh/hdf5";
      extracts["e8/params/fields"].append().set("GPhi1");
      extracts["e8/params/fields"].append().set("GPhi2");
      extracts["e8/params/fields"].append().set("GPhi3");
    
} // defineActions() call end here


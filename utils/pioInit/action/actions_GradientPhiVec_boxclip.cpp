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


void parIO::defineActions_GradientPhiVec_boxclip(glsol &sol) {

      /* this is fuction only for relay GPhi extact, not for insitu */

      conduit::Node &add_act7 = actions.append();
      add_act7["action"] = "add_pipelines";
      conduit::Node &pipelines5 = add_act7["pipelines"];
      pipelines5["pl17/f1/type"] = "clip";
      conduit::Node &clip_params = pipelines5["pl17/f1/params"];
      clip_params["topology"] = "topo";
      clip_params["box/min/x"] = sol.config.GPhi_boxclip1_min_x;
      clip_params["box/min/y"] = sol.config.GPhi_boxclip1_min_y;
      clip_params["box/min/z"] = sol.config.GPhi_boxclip1_min_z;
      clip_params["box/max/x"] = sol.config.GPhi_boxclip1_max_x;
      clip_params["box/max/y"] = sol.config.GPhi_boxclip1_max_y;
      clip_params["box/max/z"] = sol.config.GPhi_boxclip1_max_z;
  
      conduit::Node &add_act13 = actions.append();
      add_act13["action"] = "add_extracts";
      conduit::Node &extracts = add_act13["extracts"];
      extracts["e9/type"] = "relay";
      extracts["e9/pipeline"]  = "pl17";
      extracts["e9/params/path"] = "pio_Vec/dyGiLa-sim-GPhiVector-boxclip_t-%09d";
      extracts["e9/params/protocol"] = "blueprint/mesh/hdf5";
      extracts["e9/params/fields"].append().set("GPhi1");
      extracts["e9/params/fields"].append().set("GPhi2");
      extracts["e9/params/fields"].append().set("GPhi3");
    
} // defineActions() call end here


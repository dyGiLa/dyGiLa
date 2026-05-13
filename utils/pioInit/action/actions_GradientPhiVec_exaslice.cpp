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


void parIO::defineActions_GradientPhiVec_exaslice(glsol &sol) {

      /* this is fuction only for relay GPhi extact, not for insitu */

      conduit::Node &add_act5 = actions.append();
      add_act5["action"] = "add_pipelines";
      conduit::Node &pipelines3 = add_act5["pipelines"];
      pipelines3["pl15/f1/type"] = "exaslice";
      conduit::Node &slice_params = pipelines3["pl15/f1/params"];
      slice_params["point/x"] = sol.config.GPhi_slice1_point_x;
      slice_params["point/y"] = sol.config.GPhi_slice1_point_y;
      slice_params["point/z"] = sol.config.GPhi_slice1_point_z;
      slice_params["normal/x"] = sol.config.GPhi_slice1_norm_x;
      slice_params["normal/y"] = sol.config.GPhi_slice1_norm_y;
      slice_params["normal/z"] = sol.config.GPhi_slice1_norm_z;
  
      conduit::Node &add_act11 = actions.append();
      add_act11["action"] = "add_extracts";
      conduit::Node &extracts = add_act11["extracts"];
      extracts["e7/type"] = "relay";
      extracts["e7/pipeline"]  = "pl15";
      extracts["e7/params/path"] = "pio_Vec/dyGiLa-sim-GPhiVector-exaslice_t-%09d";
      extracts["e7/params/protocol"] = "blueprint/mesh/hdf5";
      extracts["e7/params/fields"].append().set("GPhi1");
      extracts["e7/params/fields"].append().set("GPhi2");
      extracts["e7/params/fields"].append().set("GPhi3");
    
} // defineActions() call end here


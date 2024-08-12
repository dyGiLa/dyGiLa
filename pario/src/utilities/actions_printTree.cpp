#define _USE_MATH_DEFINES
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


void parIO::defineActions_printTree() {

    // conduit::Node &add_act = actions.append();
    
    // add_act["action"] = "add_scenes";
    // conduit::Node &scenes = add_act["scenes"];

    // /* >>>>>>>>>>>>>>   1st scene    <<<<<<<<<<<<< */
    // scenes["s1/plots/p1/type"] = "pseudocolor";
    // scenes["s1/plots/p1/field"] = "gapAOrdered";

    // // color map clamping. min_value will be set to 0.0 if initialCondtion is 2 i.e., normal_phase
    // scenes["s1/plots/p1/min_value"]
    //   = (sol.config.initialCondition == 2 || sol.config.initialCondition == 0)
    //      ? 0.0 :  (matep.gap_A_td(sol.config.Inip, sol.config.IniT)
    // 		   * (1. + sol.config.clamp_bias_gapMin));

    // scenes["s1/plots/p1/max_value"]
    //   = matep.gap_B_td(sol.config.Inip, sol.config.IniT)
    //     * (1. + sol.config.clamp_bias_gapMax);
    
    // scenes["s1/renders/r1/image_prefix"] = "gapA_t-%04d";
    // scenes["s1/renders/r1/camera/azimuth"] = sol.config.camera_azi/*35.0*/;
    // scenes["s1/renders/r1/camera/elevation"] = sol.config.camera_ele/*30.0*/;
  
    /* >>>>>>>>>>>>>> pipleline Node <<<<<<<<<<<<< */
    
    // conduit::Node &add_act2 = actions.append();
    // add_act2["action"] = "add_pipelines";
    // conduit::Node &pipelines = add_act2["pipelines"];

    // /* >>>>>>>>>>>>>> pipleline clip <<<<<<<<<<<<< */
    
    // if (config.do_gapA_clip == 1){    
    //  pipelines["pl1/f1/type"] = "clip";
    //  conduit::Node &clip_params = pipelines["pl1/f1/params"];

    //  clip_params["topology"] = "topo";
    //  clip_params["plane/point/x"] = 40.;
    //  clip_params["plane/point/y"] = 32.;
    //  clip_params["plane/point/z"] = 32.;
    //  clip_params["plane/normal/x"] = 0.;
    //  clip_params["plane/normal/y"] = 0.;
    //  clip_params["plane/normal/z"] = -1.;

    //  scenes["s2/plots/p1/type"] = "pseudocolor";
    //  scenes["s2/plots/p1/pipeline"] = "pl1";
    //  scenes["s2/plots/p1/field"] = "gapAOrdered";

    //  scenes["s2/plots/p1/min_value"] = (config.initialCondition == 2 || config.initialCondition == 0)
    //                                    ? 0.0 :  matep.gap_A_td(config.p, config.T) * (1. + config.clamp_bias_gapMin);
    //  scenes["s2/plots/p1/max_value"] = matep.gap_B_td(config.p, config.T) + config.clamp_bias_gapMax;
    
    //  scenes["s2/renders/r1/image_prefix"] = "gapA-clip_t-%04d";
    //  scenes["s2/renders/r1/camera/azimuth"] = config.camera_azi/*35.0*/;
    //  scenes["s2/renders/r1/camera/elevation"] = config.camera_ele/*30.0*/;
    // }

    // /* >>>>>>>>>>> pipleline isosurfece <<<<<<<<<<<<< */

    // if (config.do_gapA_isosurface == 1){
    //  pipelines["pl2/f1/type"] = "contour";

    //  conduit::Node &contour_params = pipelines["pl2/f1/params"];
    //  contour_params["field"] = "gapAOrdered";
    // //double iso_vals[21] = {3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4.0, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 4.5, 4.55, 4.6, 4.65, 4.7};
    // //double iso_vals[2] = {3.45, 3.55};
    // //double iso_vals[21] = {3.20, 3.21, 3.22, 3.23, 3.24, 3.25, 3.26, 3.27, 3.28, 3.29, 3.30, 3.31, 3.32, 3.33, 3.34, 3.35, 3.36, 3.37, 3.38, 3.39, 3.40};
    //  double iso_vals[7] = {2.8, 2.90, 3.10, 3.20, 3.30, 3.40, 3.50}; 
    // //double iso_vals[11] = {2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0};
    // //double iso_vals[12] = {2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9};

    //  contour_params["iso_values"].set(iso_vals, 7);

    //  scenes["s3/plots/p1/type"] = "pseudocolor";
    //  scenes["s3/plots/p1/pipeline"] = "pl2";
    //  scenes["s3/plots/p1/field"] = "gapAOrdered";
    //  scenes["s3/renders/r1/image_prefix"] = "gapA-iso_t-%04d";

    //  double box_bounds[6] = {0.0, config.lx * config.dx, 0.0, config.ly * config.dx, 0.0, config.lz * config.dx};
    //  scenes["s3/renders/r1/dataset_bounds"].set(box_bounds,6);
    //  scenes["s3/renders/r1/render_bg"] = "true";
    
    //  scenes["s3/renders/r1/camera/azimuth"] = config.camera_azi/*35.0*/;
    //  scenes["s3/renders/r1/camera/elevation"] = config.camera_ele/*30.0*/;
    // }

    // /* >>>>>>>>>>>>>>   3slice  <<<<<<<<<<<<<<<< */
    // if (config.do_gapA_3slice == 1){
    //  pipelines["pl3/f1/type"] = "3slice";

    //  conduit::Node &slice3_params = pipelines["pl3/f1/params"];
    //  slice3_params["x_offset"] = -0.531f;
    //  slice3_params["y_offset"] = -0.6875f;
    //  slice3_params["z_offset"] = -0.21f;

    //  scenes["s4/plots/p1/type"] = "pseudocolor";
    //  scenes["s4/plots/p1/pipeline"] = "pl3";
    //  scenes["s4/plots/p1/field"] = "gapAOrdered";
    //  scenes["s4/renders/r1/image_prefix"] = "gapA-3slice_t-%04d";

    //  scenes["s4/plots/p1/min_value"] = 0.;
    //  scenes["s4/plots/p1/max_value"] = matep.gap_B_td(config.p, config.T) + config.clamp_bias_gapMax;
  
    //  scenes["s4/renders/r1/camera/azimuth"] = config.camera_azi/*35.0*/;
    //  scenes["s4/renders/r1/camera/elevation"] = config.camera_ele/*30.0*/;
    // }
    // /* >>>>>>>>>>>>>>   slice     <<<<<<<<<<<<< */
    // if (config.do_fe_slice == 1){
    //   pipelines["pl4/f1/type"] = "slice";
    //   // filter knobs
    //   conduit::Node &slice_params = pipelines["pl4/f1/params"];
    //   slice_params["point/x"] = 100.f;
    //   slice_params["point/y"] = 50.f;
    //   slice_params["point/z"] = 60.f;

    //   slice_params["normal/x"] = 0.f;
    //   slice_params["normal/y"] = 0.f;
    //   slice_params["normal/z"] = 1.f;

    //   scenes["s6/plots/p1/type"] = "pseudocolor";
    //   scenes["s6/plots/p1/pipeline"] = "pl4";
    //   scenes["s6/plots/p1/field"] = "feDensityOrdered";
    //   scenes["s6/renders/r1/image_prefix"] = "FED-slice_t-%04d";

    //   scenes["s6/plots/p1/min_value"] = config.clamp_fed_Min;
    //   scenes["s6/plots/p1/max_value"] = config.clamp_fed_Max;
  
    //   scenes["s6/renders/r1/camera/azimuth"] = 0.0/*35.0*/;
    //   scenes["s6/renders/r1/camera/elevation"] = 0.0/*30.0*/;      
    // }

    // if (config.do_gapA_slice == 1){
    //   pipelines["pl5/f1/type"] = "slice";
    //   // filter knobs
    //   conduit::Node &slice_params = pipelines["pl5/f1/params"];
    //   slice_params["point/x"] = 100.f;
    //   slice_params["point/y"] = 50.f;
    //   slice_params["point/z"] = 60.f;

    //   slice_params["normal/x"] = 0.f;
    //   slice_params["normal/y"] = 0.f;
    //   slice_params["normal/z"] = 1.f;

    //   scenes["s7/plots/p1/type"] = "pseudocolor";
    //   scenes["s7/plots/p1/pipeline"] = "pl5";
    //   scenes["s7/plots/p1/field"] = "gapAOrdered";
    //   scenes["s7/renders/r1/image_prefix"] = "gapA-slice_t-%04d";

    //   scenes["s7/plots/p1/min_value"] = 0.;
    //   scenes["s7/plots/p1/max_value"] = 4.5;
  
    //   scenes["s7/renders/r1/camera/azimuth"] = 0.0/*35.0*/;
    //   scenes["s7/renders/r1/camera/elevation"] = 0.0/*30.0*/;      
    // }
    
    // /* >>>>>>>>>>>>>>   2nd scene    <<<<<<<<<<<<< */    
    // scenes["s5/plots/p1/type"] = "pseudocolor";
    // scenes["s5/plots/p1/field"] = "feDensityOrdered";

    // scenes["s5/plots/p1/min_value"] = 0.0;					  
    // scenes["s5/plots/p1/max_value"] = 5.0;
    
    // scenes["s5/renders/r1/image_prefix"] = "feDensity_t-%04d";
    // scenes["s5/renders/r1/camera/azimuth"] = config.camera_azi/*35.0*/;
    // scenes["s5/renders/r1/camera/elevation"] = config.camera_ele/*30.0*/;
   
    // /* >>>>>>>>>>>>> extract hdf5 <<<<<<<<<<<<<< */

    // if (config.hdf5_A_matrix_output == 1){
    //  conduit::Node &add_act3 = actions.append();
    //  add_act3["action"] = "add_extracts";

    //  conduit::Node &extracts = add_act3["extracts"];
    //  extracts["e1/type"] = "relay";
    //  extracts["e1/params/path"] = "sim-data";
    //  extracts["e1/params/protocol"] = "blueprint/mesh/hdf5";

    //  extracts["e1/params/fields"].append().set("gapAOrdered");
    //  extracts["e1/params/fields"].append().set("feDensityOrdered");
    
    //  extracts["e1/params/fields"].append().set("u11Ordered");
    //  extracts["e1/params/fields"].append().set("u12Ordered");
    //  extracts["e1/params/fields"].append().set("u13Ordered");
    //  extracts["e1/params/fields"].append().set("u21Ordered");
    //  extracts["e1/params/fields"].append().set("u22Ordered");
    //  extracts["e1/params/fields"].append().set("u23Ordered");
    //  extracts["e1/params/fields"].append().set("u31Ordered");
    //  extracts["e1/params/fields"].append().set("u32Ordered");
    //  extracts["e1/params/fields"].append().set("u33Ordered");

    //  extracts["e1/params/fields"].append().set("v11Ordered");
    //  extracts["e1/params/fields"].append().set("v12Ordered");
    //  extracts["e1/params/fields"].append().set("v13Ordered");
    //  extracts["e1/params/fields"].append().set("v21Ordered");
    //  extracts["e1/params/fields"].append().set("v22Ordered");
    //  extracts["e1/params/fields"].append().set("v23Ordered");
    //  extracts["e1/params/fields"].append().set("v31Ordered");
    //  extracts["e1/params/fields"].append().set("v32Ordered");
    //  extracts["e1/params/fields"].append().set("v33Ordered");
    // }

    // if (config.hdf5_trA_output == 1){
    //  conduit::Node &add_act4 = actions.append();
    //  add_act4["action"] = "add_extracts";

    //  conduit::Node &extracts = add_act4["extracts"];
    //  extracts["e1/type"] = "relay";
    //  extracts["e1/params/path"] = "sim-data";
    //  extracts["e1/params/protocol"] = "blueprint/mesh/hdf5";

    //  extracts["e1/params/fields"].append().set("gapAOrdered");
    //  extracts["e1/params/fields"].append().set("feDensityOrdered");
    //  extracts["e1/params/fields"].append().set("trA_reOrdered");
    //  extracts["e1/params/fields"].append().set("trA_imOrdered");     
 
    // }

    // if (config.hdf5_eigvA_output == 1){
    //   conduit::Node &add_act5 = actions.append();
    //   add_act5["action"] = "add_extracts";

    //   conduit::Node &extracts = add_act5["extracts"];
    //   extracts["e1/type"] = "relay";
    //   extracts["e1/params/path"] = "sim-data";
    //   extracts["e1/params/protocol"] = "blueprint/mesh/hdf5";

    //   extracts["e1/params/fields"].append().set("gapAOrdered");
    //   extracts["e1/params/fields"].append().set("feDensityOrdered");
    //   extracts["e1/params/fields"].append().set("eigAv1Ordered");
    //   extracts["e1/params/fields"].append().set("eigAv2Ordered");
    //   extracts["e1/params/fields"].append().set("eigAv3Ordered");           
    // }

    // if (config.hdf5_mass_current_output == 1){
    //   conduit::Node &add_act6 = actions.append();
    //   add_act6["action"] = "add_extracts";

    //   conduit::Node &extracts = add_act6["extracts"];
    //   extracts["e1/type"] = "relay";
    //   extracts["e1/params/path"] = "sim-data";
    //   extracts["e1/params/protocol"] = "blueprint/mesh/hdf5";

    //   extracts["e1/params/fields"].append().set("gapAOrdered");
    //   extracts["e1/params/fields"].append().set("feDensityOrdered");
    //   extracts["e1/params/fields"].append().set("jm1Ordered");
    //   extracts["e1/params/fields"].append().set("jm2Ordered");
    //   extracts["e1/params/fields"].append().set("jm3Ordered");
    //   extracts["e1/params/fields"].append().set("phaseExpModulusO");
    //   extracts["e1/params/fields"].append().set("phaseExpAngleO");
    //   extracts["e1/params/fields"].append().set("phaseExp2ReO");
    //   extracts["e1/params/fields"].append().set("phaseExp2ImO");                       
    // }

    // if (config.hdf5_spin_current_output == 1){
    //   conduit::Node &add_act7 = actions.append();
    //   add_act7["action"] = "add_extracts";

    //   conduit::Node &extracts = add_act7["extracts"];
    //   extracts["e1/type"] = "relay";
    //   extracts["e1/params/path"] = "sim-data";
    //   extracts["e1/params/protocol"] = "blueprint/mesh/hdf5";

    //   extracts["e1/params/fields"].append().set("gapAOrdered");
    //   extracts["e1/params/fields"].append().set("feDensityOrdered");
    //   extracts["e1/params/fields"].append().set("js11O");
    //   extracts["e1/params/fields"].append().set("js21O");
    //   extracts["e1/params/fields"].append().set("js31O");
    //   extracts["e1/params/fields"].append().set("js12O");
    //   extracts["e1/params/fields"].append().set("js22O");
    //   extracts["e1/params/fields"].append().set("js32O");
    //   extracts["e1/params/fields"].append().set("js13O");
    //   extracts["e1/params/fields"].append().set("js23O");
    //   extracts["e1/params/fields"].append().set("js33O");      
    // }
    
    
    /* >>>>>>>>>>>>> ???????????? <<<<<<<<<<<<<< */
     
    // print our full actions tree
    hila::out0 << actions.to_yaml()
	       << '\n'
	       << std::endl;
    
} // defineActions() call end here


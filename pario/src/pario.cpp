#define _USE_MATH_DEFINES
#define USE_ASCENT 
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
//#include "matep.hpp"
#include "pario.hpp"

/*-----------------------------------------------------------------------*/
/***      include Ascent & Canduit for in situ rank rendering        *****/
/*-----------------------------------------------------------------------*/
#if defined USE_ASCENT

#include "ascent.hpp"
#include "conduit_blueprint.hpp"

#endif
/*-----------------------------------------------------------------------*/
/***               include Ascent & Canduit end here                 *****/
/*-----------------------------------------------------------------------*/

void he3sim::insitu_createMesh() {

    // Create a 3D mesh defined on a uniform grid of points
   
    // conduit::Node mesh;
    mesh["state/time"].set_external(&t);
    // mesh["state/cycle"].set_external(&step);
#if defined USE_MPI
    mesh["state/domain_id"] = lattice.mynode.rank;
#endif
    mesh["state/software"] = "He3Sim";
    mesh["state/title"] = "Bulk TDGL equation simulator";
    mesh["state/info"] = "In Situ rendering of order parameter field from He3Sim";

    // create the coordinate set
    mesh["coordsets/coords/type"] = "uniform";
    mesh["coordsets/coords/dims/i"] = lattice.mynode.size[0] + 2;
    mesh["coordsets/coords/dims/j"] = lattice.mynode.size[1] + 2;
    mesh["coordsets/coords/dims/k"] = lattice.mynode.size[2] + 2;

    // add origin and spacing to the coordset (optional)
    mesh["coordsets/coords/origin/x"] = ((lattice.mynode.min[0] - 1) * config.dx);
    mesh["coordsets/coords/origin/y"] = ((lattice.mynode.min[1] - 1) * config.dx);
    mesh["coordsets/coords/origin/z"] = ((lattice.mynode.min[2] - 1) * config.dx);
    mesh["coordsets/coords/spacing/dx"] = config.dx;
    mesh["coordsets/coords/spacing/dy"] = config.dx;
    mesh["coordsets/coords/spacing/dz"] = config.dx;
    
    // add the topology
    // this case is simple b/c it's implicitly derived from the coordinate set
    mesh["topologies/topo/type"] = "uniform";
    // reference the coordinate set by name
    mesh["topologies/topo/coordset"] = "coords";

    // create an vertex associated field named gapAOrdered
    mesh["fields/gapAOrdered/association"] = "vertex";
    mesh["fields/gapAOrdered/topology"] = "topo";
    mesh["fields/gapAOrdered/values"].set_external(gapAOrdered.data(), latticeVolumeWithGhost);

    // create an vertex associated field named feDensityOrdered
    mesh["fields/feDensityOrdered/association"] = "vertex";
    mesh["fields/feDensityOrdered/topology"] = "topo";
    mesh["fields/feDensityOrdered/values"].set_external(feDensityOrdered.data(), latticeVolumeWithGhost);

    // create an vertex associated field named trAOrdered
    if (config.hdf5_trA_output == 1){
      mesh["fields/trA_reOrdered/association"] = "vertex";
      mesh["fields/trA_reOrdered/topology"] = "topo";
      mesh["fields/trA_reOrdered/values"].set_external(trA_reOrdered.data(), latticeVolumeWithGhost);

      mesh["fields/trA_imOrdered/association"] = "vertex";
      mesh["fields/trA_imOrdered/topology"] = "topo";
      mesh["fields/trA_imOrdered/values"].set_external(trA_imOrdered.data(), latticeVolumeWithGhost);      
    }

    if (config.hdf5_eigvA_output == 1){
      mesh["fields/eigAv1Ordered/association"] = "vertex";
      mesh["fields/eigAv1Ordered/topology"] = "topo";
      mesh["fields/eigAv1Ordered/values"].set_external(eigAv1Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/eigAv2Ordered/association"] = "vertex";
      mesh["fields/eigAv2Ordered/topology"] = "topo";
      mesh["fields/eigAv2Ordered/values"].set_external(eigAv2Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/eigAv3Ordered/association"] = "vertex";
      mesh["fields/eigAv3Ordered/topology"] = "topo";
      mesh["fields/eigAv3Ordered/values"].set_external(eigAv3Ordered.data(), latticeVolumeWithGhost);      
    }

    if (config.hdf5_mass_current_output == 1){
      mesh["fields/jm1Ordered/association"] = "vertex";
      mesh["fields/jm1Ordered/topology"] = "topo";
      mesh["fields/jm1Ordered/values"].set_external(jm1Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/jm2Ordered/association"] = "vertex";
      mesh["fields/jm2Ordered/topology"] = "topo";
      mesh["fields/jm2Ordered/values"].set_external(jm2Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/jm3Ordered/association"] = "vertex";
      mesh["fields/jm3Ordered/topology"] = "topo";
      mesh["fields/jm3Ordered/values"].set_external(jm3Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/phaseExpModulusO/association"] = "vertex";
      mesh["fields/phaseExpModulusO/topology"] = "topo";
      mesh["fields/phaseExpModulusO/values"].set_external(phaseExpModulusO.data(), latticeVolumeWithGhost);

      mesh["fields/phaseExpAngleO/association"] = "vertex";
      mesh["fields/phaseExpAngleO/topology"] = "topo";
      mesh["fields/phaseExpAngleO/values"].set_external(phaseExpAngleO.data(), latticeVolumeWithGhost);

      mesh["fields/phaseExp2ReO/association"] = "vertex";
      mesh["fields/phaseExp2ReO/topology"] = "topo";
      mesh["fields/phaseExp2ReO/values"].set_external(phaseExp2ReO.data(), latticeVolumeWithGhost);

      mesh["fields/phaseExp2ImO/association"] = "vertex";
      mesh["fields/phaseExp2ImO/topology"] = "topo";
      mesh["fields/phaseExp2ImO/values"].set_external(phaseExp2ImO.data(), latticeVolumeWithGhost);      
    }    

    if (config.hdf5_spin_current_output == 1){
      mesh["fields/js11O/association"] = "vertex";
      mesh["fields/js11O/topology"] = "topo";
      mesh["fields/js11O/values"].set_external(js11O.data(), latticeVolumeWithGhost);

      mesh["fields/js21O/association"] = "vertex";
      mesh["fields/js21O/topology"] = "topo";
      mesh["fields/js21O/values"].set_external(js21O.data(), latticeVolumeWithGhost);

      mesh["fields/js31O/association"] = "vertex";
      mesh["fields/js31O/topology"] = "topo";
      mesh["fields/js31O/values"].set_external(js31O.data(), latticeVolumeWithGhost);

      mesh["fields/js12O/association"] = "vertex";
      mesh["fields/js12O/topology"] = "topo";
      mesh["fields/js12O/values"].set_external(js12O.data(), latticeVolumeWithGhost);

      mesh["fields/js22O/association"] = "vertex";
      mesh["fields/js22O/topology"] = "topo";
      mesh["fields/js22O/values"].set_external(js22O.data(), latticeVolumeWithGhost);

      mesh["fields/js32O/association"] = "vertex";
      mesh["fields/js32O/topology"] = "topo";
      mesh["fields/js32O/values"].set_external(js32O.data(), latticeVolumeWithGhost);

      mesh["fields/js13O/association"] = "vertex";
      mesh["fields/js13O/topology"] = "topo";
      mesh["fields/js13O/values"].set_external(js13O.data(), latticeVolumeWithGhost);

      mesh["fields/js23O/association"] = "vertex";
      mesh["fields/js23O/topology"] = "topo";
      mesh["fields/js23O/values"].set_external(js23O.data(), latticeVolumeWithGhost);

      mesh["fields/js33O/association"] = "vertex";
      mesh["fields/js33O/topology"] = "topo";
      mesh["fields/js33O/values"].set_external(js33O.data(), latticeVolumeWithGhost);
   
    }    

    
    /*----------------------------------------------------------------------*/
    /*---create vertices associated field named of uxxOrdered vxxOrdered ---*/
    /*----------------------------------------------------------------------*/
    if (config.A_matrix_output == 1){
      mesh["fields/u11Ordered/association"] = "vertex";
      mesh["fields/u11Ordered/topology"] = "topo";
      mesh["fields/u11Ordered/values"].set_external(u11Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/u12Ordered/association"] = "vertex";
      mesh["fields/u12Ordered/topology"] = "topo";
      mesh["fields/u12Ordered/values"].set_external(u12Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/u13Ordered/association"] = "vertex";
      mesh["fields/u13Ordered/topology"] = "topo";
      mesh["fields/u13Ordered/values"].set_external(u13Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/u21Ordered/association"] = "vertex";
      mesh["fields/u21Ordered/topology"] = "topo";
      mesh["fields/u21Ordered/values"].set_external(u21Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/u22Ordered/association"] = "vertex";
      mesh["fields/u22Ordered/topology"] = "topo";
      mesh["fields/u22Ordered/values"].set_external(u22Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/u23Ordered/association"] = "vertex";
      mesh["fields/u23Ordered/topology"] = "topo";
      mesh["fields/u23Ordered/values"].set_external(u23Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/u31Ordered/association"] = "vertex";
      mesh["fields/u31Ordered/topology"] = "topo";
      mesh["fields/u31Ordered/values"].set_external(u31Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/u32Ordered/association"] = "vertex";
      mesh["fields/u32Ordered/topology"] = "topo";
      mesh["fields/u32Ordered/values"].set_external(u32Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/u33Ordered/association"] = "vertex";
      mesh["fields/u33Ordered/topology"] = "topo";
      mesh["fields/u33Ordered/values"].set_external(u33Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/v11Ordered/association"] = "vertex";
      mesh["fields/v11Ordered/topology"] = "topo";
      mesh["fields/v11Ordered/values"].set_external(v11Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/v12Ordered/association"] = "vertex";
      mesh["fields/v12Ordered/topology"] = "topo";
      mesh["fields/v12Ordered/values"].set_external(v12Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/v13Ordered/association"] = "vertex";
      mesh["fields/v13Ordered/topology"] = "topo";
      mesh["fields/v13Ordered/values"].set_external(v13Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/v21Ordered/association"] = "vertex";
      mesh["fields/v21Ordered/topology"] = "topo";
      mesh["fields/v21Ordered/values"].set_external(v21Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/v22Ordered/association"] = "vertex";
      mesh["fields/v22Ordered/topology"] = "topo";
      mesh["fields/v22Ordered/values"].set_external(v22Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/v23Ordered/association"] = "vertex";
      mesh["fields/v23Ordered/topology"] = "topo";
      mesh["fields/v23Ordered/values"].set_external(v23Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/v31Ordered/association"] = "vertex";
      mesh["fields/v31Ordered/topology"] = "topo";
      mesh["fields/v31Ordered/values"].set_external(v31Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/v32Ordered/association"] = "vertex";
      mesh["fields/v32Ordered/topology"] = "topo";
      mesh["fields/v32Ordered/values"].set_external(v32Ordered.data(), latticeVolumeWithGhost);

      mesh["fields/v33Ordered/association"] = "vertex";
      mesh["fields/v33Ordered/topology"] = "topo";
      mesh["fields/v33Ordered/values"].set_external(v33Ordered.data(), latticeVolumeWithGhost);              
    }
    /*----------------------------------------------------------------------*/
    /*--- vertex associated field named of uxxOrdered vxxOrdered end here --*/
    /*----------------------------------------------------------------------*/

    
    // create an element associated field named ghostCells
    mesh["fields/ascent_ghosts/association"] = "element";
    mesh["fields/ascent_ghosts/topology"] = "topo";
    mesh["fields/ascent_ghosts/values"].set_external(ghostCellsMask, ghostVolume);

    // make sure the mesh we created conforms to the blueprint
    conduit::Node verify_info;
    if (!conduit::blueprint::mesh::verify(mesh, verify_info)) {
        hila::out0 << "Mesh Verify failed!\n";
        hila::out0 << verify_info.to_yaml() << '\n';
    }
    else {
        hila::out0 << "Mesh verify success!\n";
    }
};

void he3sim::insitu_defineActions() {

    conduit::Node &add_act = actions.append();
    
    add_act["action"] = "add_scenes";
    conduit::Node &scenes = add_act["scenes"];

    /* >>>>>>>>>>>>>>   1st scene    <<<<<<<<<<<<< */
    scenes["s1/plots/p1/type"] = "pseudocolor";
    scenes["s1/plots/p1/field"] = "gapAOrdered";

    // color map clamping. min_value will be set to 0.0 if initialCondtion is 2 i.e., normal_phase
    scenes["s1/plots/p1/min_value"] = (config.initialCondition == 2 || config.initialCondition == 0)
                                       ? 0.0 :  matep.gap_A_td(config.p, config.T) * (1. + config.clamp_bias_gapMin);
    scenes["s1/plots/p1/max_value"] = matep.gap_B_td(config.p, config.T) * (1. + config.clamp_bias_gapMax);
    
    scenes["s1/renders/r1/image_prefix"] = "gapA_t-%04d";
    scenes["s1/renders/r1/camera/azimuth"] = config.camera_azi/*35.0*/;
    scenes["s1/renders/r1/camera/elevation"] = config.camera_ele/*30.0*/;
  
    /* >>>>>>>>>>>>>> pipleline Node <<<<<<<<<<<<< */
    
    conduit::Node &add_act2 = actions.append();
    add_act2["action"] = "add_pipelines";
    conduit::Node &pipelines = add_act2["pipelines"];

    /* >>>>>>>>>>>>>> pipleline clip <<<<<<<<<<<<< */
    
    if (config.do_gapA_clip == 1){    
     pipelines["pl1/f1/type"] = "clip";
     conduit::Node &clip_params = pipelines["pl1/f1/params"];

     clip_params["topology"] = "topo";
     clip_params["plane/point/x"] = 40.;
     clip_params["plane/point/y"] = 32.;
     clip_params["plane/point/z"] = 32.;
     clip_params["plane/normal/x"] = 0.;
     clip_params["plane/normal/y"] = 0.;
     clip_params["plane/normal/z"] = -1.;

     scenes["s2/plots/p1/type"] = "pseudocolor";
     scenes["s2/plots/p1/pipeline"] = "pl1";
     scenes["s2/plots/p1/field"] = "gapAOrdered";

     scenes["s2/plots/p1/min_value"] = (config.initialCondition == 2 || config.initialCondition == 0)
                                       ? 0.0 :  matep.gap_A_td(config.p, config.T) * (1. + config.clamp_bias_gapMin);
     scenes["s2/plots/p1/max_value"] = matep.gap_B_td(config.p, config.T) + config.clamp_bias_gapMax;
    
     scenes["s2/renders/r1/image_prefix"] = "gapA-clip_t-%04d";
     scenes["s2/renders/r1/camera/azimuth"] = config.camera_azi/*35.0*/;
     scenes["s2/renders/r1/camera/elevation"] = config.camera_ele/*30.0*/;
    }

    /* >>>>>>>>>>> pipleline isosurfece <<<<<<<<<<<<< */

    if (config.do_gapA_isosurface == 1){
     pipelines["pl2/f1/type"] = "contour";

     conduit::Node &contour_params = pipelines["pl2/f1/params"];
     contour_params["field"] = "gapAOrdered";
    //double iso_vals[21] = {3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4.0, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 4.5, 4.55, 4.6, 4.65, 4.7};
    //double iso_vals[2] = {3.45, 3.55};
    //double iso_vals[21] = {3.20, 3.21, 3.22, 3.23, 3.24, 3.25, 3.26, 3.27, 3.28, 3.29, 3.30, 3.31, 3.32, 3.33, 3.34, 3.35, 3.36, 3.37, 3.38, 3.39, 3.40};
     double iso_vals[7] = {2.8, 2.90, 3.10, 3.20, 3.30, 3.40, 3.50}; 
    //double iso_vals[11] = {2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0};
    //double iso_vals[12] = {2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9};

     contour_params["iso_values"].set(iso_vals, 7);

     scenes["s3/plots/p1/type"] = "pseudocolor";
     scenes["s3/plots/p1/pipeline"] = "pl2";
     scenes["s3/plots/p1/field"] = "gapAOrdered";
     scenes["s3/renders/r1/image_prefix"] = "gapA-iso_t-%04d";

     double box_bounds[6] = {0.0, config.lx * config.dx, 0.0, config.ly * config.dx, 0.0, config.lz * config.dx};
     scenes["s3/renders/r1/dataset_bounds"].set(box_bounds,6);
     scenes["s3/renders/r1/render_bg"] = "true";
    
     scenes["s3/renders/r1/camera/azimuth"] = config.camera_azi/*35.0*/;
     scenes["s3/renders/r1/camera/elevation"] = config.camera_ele/*30.0*/;
    }

    /* >>>>>>>>>>>>>>   3slice  <<<<<<<<<<<<<<<< */
    if (config.do_gapA_3slice == 1){
     pipelines["pl3/f1/type"] = "3slice";

     conduit::Node &slice3_params = pipelines["pl3/f1/params"];
     slice3_params["x_offset"] = -0.531f;
     slice3_params["y_offset"] = -0.6875f;
     slice3_params["z_offset"] = -0.21f;

     scenes["s4/plots/p1/type"] = "pseudocolor";
     scenes["s4/plots/p1/pipeline"] = "pl3";
     scenes["s4/plots/p1/field"] = "gapAOrdered";
     scenes["s4/renders/r1/image_prefix"] = "gapA-3slice_t-%04d";

     scenes["s4/plots/p1/min_value"] = 0.;
     scenes["s4/plots/p1/max_value"] = matep.gap_B_td(config.p, config.T) + config.clamp_bias_gapMax;
  
     scenes["s4/renders/r1/camera/azimuth"] = config.camera_azi/*35.0*/;
     scenes["s4/renders/r1/camera/elevation"] = config.camera_ele/*30.0*/;
    }
    /* >>>>>>>>>>>>>>   slice     <<<<<<<<<<<<< */
    if (config.do_fe_slice == 1){
      pipelines["pl4/f1/type"] = "slice";
      // filter knobs
      conduit::Node &slice_params = pipelines["pl4/f1/params"];
      slice_params["point/x"] = 100.f;
      slice_params["point/y"] = 50.f;
      slice_params["point/z"] = 60.f;

      slice_params["normal/x"] = 0.f;
      slice_params["normal/y"] = 0.f;
      slice_params["normal/z"] = 1.f;

      scenes["s6/plots/p1/type"] = "pseudocolor";
      scenes["s6/plots/p1/pipeline"] = "pl4";
      scenes["s6/plots/p1/field"] = "feDensityOrdered";
      scenes["s6/renders/r1/image_prefix"] = "FED-slice_t-%04d";

      scenes["s6/plots/p1/min_value"] = config.clamp_fed_Min;
      scenes["s6/plots/p1/max_value"] = config.clamp_fed_Max;
  
      scenes["s6/renders/r1/camera/azimuth"] = 0.0/*35.0*/;
      scenes["s6/renders/r1/camera/elevation"] = 0.0/*30.0*/;      
    }

    if (config.do_gapA_slice == 1){
      pipelines["pl5/f1/type"] = "slice";
      // filter knobs
      conduit::Node &slice_params = pipelines["pl5/f1/params"];
      slice_params["point/x"] = 100.f;
      slice_params["point/y"] = 50.f;
      slice_params["point/z"] = 60.f;

      slice_params["normal/x"] = 0.f;
      slice_params["normal/y"] = 0.f;
      slice_params["normal/z"] = 1.f;

      scenes["s7/plots/p1/type"] = "pseudocolor";
      scenes["s7/plots/p1/pipeline"] = "pl5";
      scenes["s7/plots/p1/field"] = "gapAOrdered";
      scenes["s7/renders/r1/image_prefix"] = "gapA-slice_t-%04d";

      scenes["s7/plots/p1/min_value"] = 0.;
      scenes["s7/plots/p1/max_value"] = 4.5;
  
      scenes["s7/renders/r1/camera/azimuth"] = 0.0/*35.0*/;
      scenes["s7/renders/r1/camera/elevation"] = 0.0/*30.0*/;      
    }
    
    /* >>>>>>>>>>>>>>   2nd scene    <<<<<<<<<<<<< */    
    scenes["s5/plots/p1/type"] = "pseudocolor";
    scenes["s5/plots/p1/field"] = "feDensityOrdered";

    scenes["s5/plots/p1/min_value"] = 0.0;					  
    scenes["s5/plots/p1/max_value"] = 5.0;
    
    scenes["s5/renders/r1/image_prefix"] = "feDensity_t-%04d";
    scenes["s5/renders/r1/camera/azimuth"] = config.camera_azi/*35.0*/;
    scenes["s5/renders/r1/camera/elevation"] = config.camera_ele/*30.0*/;
   
    /* >>>>>>>>>>>>> extract hdf5 <<<<<<<<<<<<<< */

    if (config.hdf5_A_matrix_output == 1){
     conduit::Node &add_act3 = actions.append();
     add_act3["action"] = "add_extracts";

     conduit::Node &extracts = add_act3["extracts"];
     extracts["e1/type"] = "relay";
     extracts["e1/params/path"] = "sim-data";
     extracts["e1/params/protocol"] = "blueprint/mesh/hdf5";

     extracts["e1/params/fields"].append().set("gapAOrdered");
     extracts["e1/params/fields"].append().set("feDensityOrdered");
    
     extracts["e1/params/fields"].append().set("u11Ordered");
     extracts["e1/params/fields"].append().set("u12Ordered");
     extracts["e1/params/fields"].append().set("u13Ordered");
     extracts["e1/params/fields"].append().set("u21Ordered");
     extracts["e1/params/fields"].append().set("u22Ordered");
     extracts["e1/params/fields"].append().set("u23Ordered");
     extracts["e1/params/fields"].append().set("u31Ordered");
     extracts["e1/params/fields"].append().set("u32Ordered");
     extracts["e1/params/fields"].append().set("u33Ordered");

     extracts["e1/params/fields"].append().set("v11Ordered");
     extracts["e1/params/fields"].append().set("v12Ordered");
     extracts["e1/params/fields"].append().set("v13Ordered");
     extracts["e1/params/fields"].append().set("v21Ordered");
     extracts["e1/params/fields"].append().set("v22Ordered");
     extracts["e1/params/fields"].append().set("v23Ordered");
     extracts["e1/params/fields"].append().set("v31Ordered");
     extracts["e1/params/fields"].append().set("v32Ordered");
     extracts["e1/params/fields"].append().set("v33Ordered");
    }

    if (config.hdf5_trA_output == 1){
     conduit::Node &add_act4 = actions.append();
     add_act4["action"] = "add_extracts";

     conduit::Node &extracts = add_act4["extracts"];
     extracts["e1/type"] = "relay";
     extracts["e1/params/path"] = "sim-data";
     extracts["e1/params/protocol"] = "blueprint/mesh/hdf5";

     extracts["e1/params/fields"].append().set("gapAOrdered");
     extracts["e1/params/fields"].append().set("feDensityOrdered");
     extracts["e1/params/fields"].append().set("trA_reOrdered");
     extracts["e1/params/fields"].append().set("trA_imOrdered");     
 
    }

    if (config.hdf5_eigvA_output == 1){
      conduit::Node &add_act5 = actions.append();
      add_act5["action"] = "add_extracts";

      conduit::Node &extracts = add_act5["extracts"];
      extracts["e1/type"] = "relay";
      extracts["e1/params/path"] = "sim-data";
      extracts["e1/params/protocol"] = "blueprint/mesh/hdf5";

      extracts["e1/params/fields"].append().set("gapAOrdered");
      extracts["e1/params/fields"].append().set("feDensityOrdered");
      extracts["e1/params/fields"].append().set("eigAv1Ordered");
      extracts["e1/params/fields"].append().set("eigAv2Ordered");
      extracts["e1/params/fields"].append().set("eigAv3Ordered");           
    }

    if (config.hdf5_mass_current_output == 1){
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
    }

    if (config.hdf5_spin_current_output == 1){
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
    }
    
    
    /* >>>>>>>>>>>>> ???????????? <<<<<<<<<<<<<< */
     
    // print our full actions tree
    hila::out0 << actions.to_yaml() << '\n' << std::endl;
}

void he3sim::insitu_hdf5xdmf(){

  const std::string fname = "rank_xmls/" + config.xmf2_fname + "_" + std::to_string(hila::myrank()) + ".xml";
  config.xml_out.open(fname, std::ios::out);

  hila::out << config.xml_out.good() << "\n";
  hila::out << config.xml_out.is_open() << "\n";
  hila::out << config.xml_out.fail() << "\n";
  
  const long dim_0 = lattice.mynode.size[0] + 2,
             dim_1 = lattice.mynode.size[1] + 2,
             dim_2 = lattice.mynode.size[2] + 2;

  unsigned int n;

    config.xml_out << "<Grid Name=\"sim-data\" Type=\"Uniform\">\n"
                    << "  <Topology name=\"topo\" TopologyType=\"3DRectMesh\" Dimensions=\""
		    << dim_2 << " " << dim_1 << " " << dim_0 << "\"" << ">" << "\n"
                    << "  </Topology>\n"
                    << "  <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n"
                    << "    <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">"
		    << "\n"
                    << "      "
                    << ((lattice.mynode.min[0] - 1) * config.dx) << " "
                    << ((lattice.mynode.min[1] - 1) * config.dx) << " "
                    << ((lattice.mynode.min[2] - 1) * config.dx) << "\n"
                    << "    </DataItem>\n"
                    << "    <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n"
                    << "      "
                    << config.dx << " " << config.dx << " " << config.dx << "\n"
                    << "    </DataItem>\n"
                    << "  </Geometry>\n"
                    << "  <Attribute Name=\"gapA\" AttributeType=\"Scalar\" Center=\"Node\">\n"
                    << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                    << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	            << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank() << ".hdf5:/fields/gapAOrdered/values"
		    << "\n"
                    << "   </DataItem>\n"
                    << "  </Attribute>\n"
                    << "  <Attribute Name=\"feDensity\" AttributeType=\"Scalar\" Center=\"Node\">\n"
                    << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                    << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	            << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank() << ".hdf5:/fields/feDensityOrdered/values"
		    << "\n"
                    << "   </DataItem>\n"
                    << "  </Attribute>"
 		    << "\n" << std::flush;
    
    if (config.hdf5_A_matrix_output == 1){
      for (n = 0; n<=8; ++n){
         config.xml_out << "  <Attribute Name=\"u"
	                 << std::to_string(n/3u + 1)
	         	 << std::to_string(n%3u + 1) << "\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/u"
		         << std::to_string(n/3u + 1)
		         << std::to_string(n%3u + 1) << "Ordered/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"
	                 << "\n"
	                 << "  <Attribute Name=\"v"
	                 << std::to_string(n/3u + 1)
		         << std::to_string(n%3u + 1) << "\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/v"
		         << std::to_string(n/3u + 1)
		         << std::to_string(n%3u + 1) << "Ordered/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"
	                 << "\n" << std::flush;
      }

    }

    if (config.hdf5_trA_output == 1) {
         config.xml_out << "  <Attribute Name=\"trA_re\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/trA_reOrdered/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"
	                 << "\n"
	                 << "  <Attribute Name=\"trA_im\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/trA_imOrdered/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"
	                 << "\n" << std::flush;
    }

    if (config.hdf5_eigvA_output == 1){
         config.xml_out << "  <Attribute Name=\"eigVal1_A\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/eigAv1Ordered/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"
	                 << "\n"
	                 << "  <Attribute Name=\"eigVal2_A\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/eigAv2Ordered/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"
			 << "\n"  
	                 << "  <Attribute Name=\"eigVal3_A\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/eigAv3Ordered/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"
			 << "\n" << std::flush;
    }

    if (config.hdf5_mass_current_output == 1){
         config.xml_out << "  <Attribute Name=\"jm_1\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/jm1Ordered/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"
	                 << "\n"
	                 << "  <Attribute Name=\"jm_2\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/jm2Ordered/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"
			 << "\n"  
	                 << "  <Attribute Name=\"jm_3\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/jm3Ordered/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"
	                 << "  <Attribute Name=\"phaseExpModulus\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/phaseExpModulusO/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"
	                 << "  <Attribute Name=\"phaseExpAngle\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/phaseExpAngleO/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"
	                 << "  <Attribute Name=\"phaseExp2ReO\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/phaseExp2ReO/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"
	                 << "  <Attribute Name=\"phaseExp2ImO\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/phaseExp2ImO/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"	   	   	   
			 << "\n" << std::flush;
    }

    if (config.hdf5_spin_current_output == 1){
          config.xml_out << "  <Attribute Name=\"js_11\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/js11O/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"
	                 << "\n"
	                 << "  <Attribute Name=\"js_21\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/js21O/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"
			 << "\n"  
	                 << "  <Attribute Name=\"js_31\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/js31O/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"
                         << "\n"	    
                         << "  <Attribute Name=\"js_12\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/js12O/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"
	                 << "\n"
	                 << "  <Attribute Name=\"js_22\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/js22O/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"
			 << "\n"  
	                 << "  <Attribute Name=\"js_32\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/js32O/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"	    
                         << "\n"
                         << "  <Attribute Name=\"js_13\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/js13O/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"
	                 << "\n"
	                 << "  <Attribute Name=\"js_23\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/js23O/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"
			 << "\n"  
	                 << "  <Attribute Name=\"js_33\""
		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                         << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                         << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		         << ".hdf5:/fields/js33O/values"
		         << "\n"
                         << "   </DataItem>\n"
                         << "  </Attribute>"	    
			 << "\n" << std::flush;
    }
     
    config.xml_out << "  <Attribute Name=\"vtkGhostType\" AttributeType=\"Scalar\" Center=\"Cell\">\n"
                    << "   <DataItem Format=\"HDF\" DataType=\"UChar\" Dimensions=\""
                    << dim_0 - 1 << " " << dim_1 - 1 << " " << dim_2 - 1 << "\"" << ">" << "\n"
 	            << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank() << ".hdf5:/fields/ascent_ghosts/values"
		    << "\n"   
                    << "   </DataItem>\n"
                    << "  </Attribute>\n"
                    << "</Grid>"
                    << "\n" << std::flush;

  config.xml_out.close();

}

void he3sim::insitu_execute() {

    /*-------------------    sqrt(Tr[A.A^+) ----------------------*/
    gapA[ALL] = real(sqrt((A[X]*A[X].dagger()).trace()));  


    /*--------------------     feDensity      --------------------*/
    real_t ebfe=fmin(matep.f_A_td(config.p, config.T), matep.f_B_td(config.p, config.T));
    feDensity[ALL] = 0;
    onsites(ALL) {
      Complex<double> a(0),b2(0),b3(0),b4(0),b5(0);
      //Complex<double> kin(0);
      Complex<double> k1(0), k2(0), k3(0);
      Complex<double> bfe(0);
      double b1 = 0;
      
      a = config.alpha * (A[X]*A[X].dagger()).trace();
      b1 = config.beta1 * ((A[X]*A[X].transpose()).trace()).squarenorm();
      b2 = config.beta2 * ((A[X]*A[X].dagger()).trace()*(A[X]*A[X].dagger()).trace());
      b3 = config.beta3 * ((A[X]*A[X].transpose()*A[X].conj()*A[X].dagger()).trace());
      b4 = config.beta4 * ((A[X]*A[X].dagger()*A[X]*A[X].dagger()).trace());
      b5 = config.beta5 * ((A[X]*A[X].dagger()*A[X].conj()*A[X].transpose()).trace());

      bfe = a + b1 + b2 + b3 + b4 + b5 - ebfe;
      //kin = (pi[X]*pi[X].dagger()).trace();
      
      foralldir(j) foralldir (k) foralldir(al){
	k1 += (A[X + k].column(j) - A[X - k].column(j)).e(al)
	      * (A[X + k].conj().column(j) - A[X - k].conj().column(j)).e(al)/(4.0*config.dx*config.dx);
	k2 += (A[X + j].column(j) - A[X - j].column(j)).e(al)
	      * (A[X + k].conj().column(k) - A[X - k].conj().column(k)).e(al)/(4.0*config.dx*config.dx);
	k3 += (A[X + k].column(j) - A[X - k].column(j)).e(al)
	      * (A[X + j].conj().column(k) - A[X - j].conj().column(k)).e(al)/(4.0*config.dx*config.dx);
      }
      // question here : how about imagnary part of k1 + k2 + k3 +bfe
      feDensity[X] = real(k1 + k2 + k3 + bfe);
    } //onsite(All) end here
    
    /*----------------     A matrix elements ---------------------*/
    if (config.A_matrix_output == 1){
     u11[ALL] = A[X].e(0,0).re; v11[ALL] = A[X].e(0,0).im;
     u12[ALL] = A[X].e(0,1).re; v12[ALL] = A[X].e(0,1).im;
     u13[ALL] = A[X].e(0,2).re; v13[ALL] = A[X].e(0,2).im;
     u21[ALL] = A[X].e(1,0).re; v21[ALL] = A[X].e(1,0).im;
     u22[ALL] = A[X].e(1,1).re; v22[ALL] = A[X].e(1,1).im;
     u23[ALL] = A[X].e(1,2).re; v23[ALL] = A[X].e(1,2).im;
     u31[ALL] = A[X].e(2,0).re; v31[ALL] = A[X].e(2,0).im;
     u32[ALL] = A[X].e(2,1).re; v32[ALL] = A[X].e(2,1).im;
     u33[ALL] = A[X].e(2,2).re; v33[ALL] = A[X].e(2,2).im;

     u11.copy_local_data_with_halo(u11Ordered); v11.copy_local_data_with_halo(v11Ordered);
     u12.copy_local_data_with_halo(u12Ordered); v12.copy_local_data_with_halo(v12Ordered);
     u13.copy_local_data_with_halo(u13Ordered); v13.copy_local_data_with_halo(v13Ordered);
     u21.copy_local_data_with_halo(u21Ordered); v21.copy_local_data_with_halo(v21Ordered);
     u22.copy_local_data_with_halo(u22Ordered); v22.copy_local_data_with_halo(v22Ordered);
     u23.copy_local_data_with_halo(u23Ordered); v23.copy_local_data_with_halo(v23Ordered);
     u31.copy_local_data_with_halo(u31Ordered); v31.copy_local_data_with_halo(v31Ordered);
     u32.copy_local_data_with_halo(u32Ordered); v32.copy_local_data_with_halo(v32Ordered);
     u33.copy_local_data_with_halo(u33Ordered); v33.copy_local_data_with_halo(v33Ordered);        
    }

    /*----------------------     trace A    ----------------------*/
    if (config.hdf5_trA_output == 1){
      trA_re[ALL] = (A[X].trace()).real();
      trA_im[ALL] = (A[X].trace()).imag();      
      trA_re.copy_local_data_with_halo(trA_reOrdered);
      trA_im.copy_local_data_with_halo(trA_imOrdered);      
    }

    /*----------------------  eigen values of A ------------------*/
    if (config.hdf5_eigvA_output == 1){
      Field<Vector<3,double>> eval;
      Field<Matrix<3,3,Complex<double>>> evec;

      onsites(ALL){
	A[X].eigen_jacobi(eval[X],evec[X]/*,hila::sort::ascending*/);
      }

      eigAv1[ALL] = eval[X].e(0); 
      eigAv2[ALL] = eval[X].e(1); 
      eigAv3[ALL] = eval[X].e(2); 

      eigAv1.copy_local_data_with_halo(eigAv1Ordered);
      eigAv2.copy_local_data_with_halo(eigAv2Ordered);
      eigAv3.copy_local_data_with_halo(eigAv3Ordered);
    }

    /*------------------ mass current components ------------------*/
    if (config.hdf5_mass_current_output == 1){
      Field<Vector<3,double>> jmX;
      Field<Complex<real_t>>  phaseExp;

      /*onsites(ALL){
        foralldir(i) foralldir(j) foralldir(al){
	  jmX[X].e(i) = ((A[X].conj().column(j)).e(al) * (A[X+i].column(j) - A[X-i].column(j)).e(al)/(2.*config.dx)
 	                + (A[X].conj().column(j)).e(al) * (A[X+j].column(i) - A[X-j].column(i)).e(al)/(2.*config.dx)
			+ (A[X].conj().column(i)).e(al) * (A[X+j].column(j) - A[X-j].column(j)).e(al)/(2.*config.dx)).imag();

        } // foralldir() calls end here

	} //onesite(ALL) call ends here*/

      onsites(ALL) {
	jmX[X] = 0;
	foralldir(i) foralldir(j) foralldir(al) {
	  jmX[X].e(i) += (A[X].e(al,j).conj() *(A[X+i].e(al,j) - A[X-i].e(al,j))
			  + A[X].e(al,j).conj() * (A[X+j].e(al,i) - A[X-j].e(al,i))
			  + A[X].e(al,i).conj() * (A[X+j].e(al,j) - A[X-j].e(al,j))).imag();
	} // foralldir end here, outermost foralldir slowest, inner run earier
	jmX[X] /= 2*config.dx;
      } // onsites(ALL) end here

      jm1[ALL] = jmX[X].e(0);
      jm2[ALL] = jmX[X].e(1);
      jm3[ALL] = jmX[X].e(2);

      /* >>>>>>>> Modulus, phase angle, phaseExp   <<<<<<< */
      phaseExp[ALL]        = ((A[X].transpose()) * A[X]).trace();
      phaseExp2Re[ALL]    = phaseExp[X].real();
      phaseExp2Im[ALL]    = phaseExp[X].imag();      
      
      phaseExpModulus[ALL] = phaseExp[X].abs();
      //phaseExpAngle[ALL]   = phaseExp[X].arg()/2.;
      phaseExpAngle[ALL]   = std::atan2(phaseExp[X].imag(), phaseExp[X].real())/2.;
      /* >>>>>>>> phase angle and modulus end here   <<<<<<*/

      jm1.copy_local_data_with_halo(jm1Ordered);
      jm2.copy_local_data_with_halo(jm2Ordered);
      jm3.copy_local_data_with_halo(jm3Ordered);

      phaseExpModulus.copy_local_data_with_halo(phaseExpModulusO);
      phaseExpAngle.copy_local_data_with_halo(phaseExpAngleO);
      phaseExp2Re.copy_local_data_with_halo(phaseExp2ReO);
      phaseExp2Im.copy_local_data_with_halo(phaseExp2ImO);            
    }

    /*------------------ spin current components ------------------*/
    if (config.hdf5_spin_current_output == 1){
      Field<Matrix<3,3,double>> jsX; // column is alpha for spin, row is i for spatial

      onsites(ALL) {
	jsX[X] = 0;
	foralldir(al) foralldir(i) foralldir(be) foralldir(ga) foralldir(j) {
          jsX[X].e(i,al) += -matep.epsilon(al,be,ga)
	                     * (A[X].e(be,i).conj() * (A[X+j].e(ga,j) - A[X-j].e(ga,j))
			        + A[X].e(be,j).conj() * (A[X+i].e(ga,j) - A[X-i].e(ga,j))
			        + A[X].e(be,j).conj() * (A[X+j].e(ga,i) - A[X-j].e(ga,i))).real();
	  
	} // foralldir end here, outermost foralldir slowest, inner run earier
	jsX[X] /= 2*config.dx;
      } // onsites(ALL) end here

      js11[ALL] = jsX[X].e(0,0);
      js21[ALL] = jsX[X].e(1,0);
      js31[ALL] = jsX[X].e(2,0);
      js12[ALL] = jsX[X].e(0,1);
      js22[ALL] = jsX[X].e(1,1);
      js32[ALL] = jsX[X].e(2,1);
      js13[ALL] = jsX[X].e(0,2);
      js23[ALL] = jsX[X].e(1,2);
      js33[ALL] = jsX[X].e(2,2);
      
      js11.copy_local_data_with_halo(js11O);
      js21.copy_local_data_with_halo(js21O);
      js31.copy_local_data_with_halo(js31O);
      
      js12.copy_local_data_with_halo(js12O);
      js22.copy_local_data_with_halo(js22O);
      js32.copy_local_data_with_halo(js32O);
      
      js13.copy_local_data_with_halo(js13O);
      js23.copy_local_data_with_halo(js23O);
      js33.copy_local_data_with_halo(js33O);

    }
    
    /*------------------------------------------------------------*/
        
    gapA.copy_local_data_with_halo(gapAOrdered);
    feDensity.copy_local_data_with_halo(feDensityOrdered);

  /* ToDo list :
   * > orbital vectors;
   * > spin vectors;
   * x mass current (done);
   * > spin currents;
   * x free energy density (done)
   * > ...
   */
    
  insitu.execute(actions);
}

void he3sim::insitu_initialize() {

    latticeVolumeWithGhost =
        (lattice.mynode.size[0] + 2) * (lattice.mynode.size[1] + 2) * (lattice.mynode.size[2] + 2);
    
    latticeVolume =
      (lattice.mynode.size[0]) * (lattice.mynode.size[1]) * (lattice.mynode.size[2]);

    gapAOrdered.reserve(latticeVolumeWithGhost);
    feDensityOrdered.reserve(latticeVolumeWithGhost);

    if (config.A_matrix_output == 1){
     u11Ordered.reserve(latticeVolumeWithGhost); v11Ordered.reserve(latticeVolumeWithGhost);
     u12Ordered.reserve(latticeVolumeWithGhost); v12Ordered.reserve(latticeVolumeWithGhost);
     u13Ordered.reserve(latticeVolumeWithGhost); v13Ordered.reserve(latticeVolumeWithGhost);
     u21Ordered.reserve(latticeVolumeWithGhost); v21Ordered.reserve(latticeVolumeWithGhost);
     u22Ordered.reserve(latticeVolumeWithGhost); v22Ordered.reserve(latticeVolumeWithGhost);
     u23Ordered.reserve(latticeVolumeWithGhost); v23Ordered.reserve(latticeVolumeWithGhost);
     u31Ordered.reserve(latticeVolumeWithGhost); v31Ordered.reserve(latticeVolumeWithGhost);
     u32Ordered.reserve(latticeVolumeWithGhost); v32Ordered.reserve(latticeVolumeWithGhost);
     u33Ordered.reserve(latticeVolumeWithGhost); v33Ordered.reserve(latticeVolumeWithGhost);    
    }

    if (config.hdf5_trA_output == 1){
      trA_reOrdered.reserve(latticeVolumeWithGhost);
      trA_imOrdered.reserve(latticeVolumeWithGhost);      
    }

    if (config.hdf5_eigvA_output == 1){
      eigAv1Ordered.reserve(latticeVolumeWithGhost);
      eigAv2Ordered.reserve(latticeVolumeWithGhost);
      eigAv3Ordered.reserve(latticeVolumeWithGhost);      
    }

    if (config.hdf5_mass_current_output == 1){
      jm1Ordered.reserve(latticeVolumeWithGhost);
      jm2Ordered.reserve(latticeVolumeWithGhost);
      jm3Ordered.reserve(latticeVolumeWithGhost);

      phaseExpModulusO.reserve(latticeVolumeWithGhost);
      phaseExpAngleO.reserve(latticeVolumeWithGhost);
      phaseExp2ReO.reserve(latticeVolumeWithGhost);
      phaseExp2ImO.reserve(latticeVolumeWithGhost);                  
    }

    if (config.hdf5_spin_current_output == 1){
      js11O.reserve(latticeVolumeWithGhost);
      js21O.reserve(latticeVolumeWithGhost);
      js31O.reserve(latticeVolumeWithGhost);

      js12O.reserve(latticeVolumeWithGhost);
      js22O.reserve(latticeVolumeWithGhost);
      js32O.reserve(latticeVolumeWithGhost);

      js13O.reserve(latticeVolumeWithGhost);
      js23O.reserve(latticeVolumeWithGhost);
      js33O.reserve(latticeVolumeWithGhost);      
    }

    
    // One more point in each direction, but cell data (Npts - 1 cells)
    auto ghostNX = lattice.mynode.size[0] + 2 - 1;
    auto ghostNY = lattice.mynode.size[1] + 2 - 1;
    auto ghostNZ = lattice.mynode.size[2] + 2 - 1;

    ghostVolume = ghostNX * ghostNY * ghostNZ;
    ghostCellsMask = (unsigned char *)memalloc(ghostVolume * sizeof(unsigned char));

    long long counter = 0;
    unsigned char Mask = 0;
    for (auto k = 0; k < ghostNZ; k++) {
        for (auto j = 0; j < ghostNY; j++) {
            for (auto i = 0; i < ghostNX; i++) {
                bool kGhostFlag = (k == 0);
                bool jGhostFlag = (j == 0);
                bool iGhostFlag = (i == 0);
                Mask = (iGhostFlag || jGhostFlag || kGhostFlag);
                ghostCellsMask[counter] = Mask;
                counter++;
            }
        }
    }

    insitu_createMesh();

    ascent_options["mpi_comm"] = MPI_Comm_c2f(lattice.mpi_comm_lat);
    ascent_options["runtime/type"] = "ascent";
#if defined CUDA
    ascent_options["runtime/vtkm/backend"] = "cuda";
    ascent_options["cuda/init"] = "false";
#endif
    ascent_options["timings"] = "false";
    
    insitu.open(ascent_options);
    insitu.publish(mesh);
    insitu_defineActions();
}

void he3sim::insitu_close() {
    insitu.close();
}

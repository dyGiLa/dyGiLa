#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>

#include "plumbing/hila.h"
#include "plumbing/fft.h"

#include "glsol.hpp"
#include "matep.hpp"
#include "pario.hpp"
#include "hsvCyclicSchemeRGBA.hpp"

#include "ascent.hpp"
#include "conduit_blueprint.hpp"


void parIO::defineActions_insitu(glsol &sol) {

    conduit::Node &add_act = actions.append();    
    add_act["action"] = "add_scenes";
    conduit::Node &scenes = add_act["scenes"];

    conduit::Node &add_act2 = actions.append();    
    add_act2["action"] = "add_scenes";
    conduit::Node &scenes2 = add_act2["scenes"];

    conduit::Node &add_act3 = actions.append();
    add_act3["action"] = "add_pipelines";
    conduit::Node &pipelines = add_act3["pipelines"];

    conduit::Node &add_act4 = actions.append();
    add_act4["action"] = "add_pipelines";
    conduit::Node &pipelines2 = add_act4["pipelines"];
    
    
    // background color vecctor
    double bg_colvec[3] = {0.8, 0.8, 0.8};
    // foreground color vector
    double fg_colvec[3] = {0., 0., 0.};
    
    /* >>>>>>>>>>>>>> pipleline gapA  clip <<<<<<<<<<<<< */
    
    if (sol.config.do_gapA_clip == 1)
      {
       pipelines["pl1/f1/type"] = "clip";
       conduit::Node &clip_params = pipelines["pl1/f1/params"];

       clip_params["topology"] = "topo";
       clip_params["plane/point/x"] = sol.config.gapA_clip1_point_x;
       clip_params["plane/point/y"] = sol.config.gapA_clip1_point_y;
       clip_params["plane/point/z"] = sol.config.gapA_clip1_point_z;
       clip_params["plane/normal/x"] = sol.config.gapA_clip1_norm_x;
       clip_params["plane/normal/y"] = sol.config.gapA_clip1_norm_y;
       clip_params["plane/normal/z"] = sol.config.gapA_clip1_norm_z;

       scenes["s2/plots/p1/type"] = "pseudocolor";
       scenes["s2/plots/p1/pipeline"] = "pl1";
       scenes["s2/plots/p1/field"] = "gapA";

       scenes["s2/plots/p1/min_value"]
	 = matep.gap_A_td(sol.config.Inip, (sol.config.Ttdb0 * matep.Tcp_mK(sol.config.Inip))) * (1. + sol.config.clamp_bias_gapMin);
    
       scenes["s2/plots/p1/max_value"]
	 = matep.gap_B_td(sol.config.Inip, (sol.config.Ttdb0 * matep.Tcp_mK(sol.config.Inip))) * (1. + sol.config.clamp_bias_gapMax);

       scenes["s2/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s2/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);    
       scenes["s2/renders/r1/image_prefix"] = "insitu/gapA-clip1/gapA-clip1_t-%09d";
       scenes["s2/renders/r1/camera/azimuth"] = sol.config.camera1_azi;
       scenes["s2/renders/r1/camera/elevation"] = sol.config.camera1_ele;
      }

    /* >>>>>>>>>>>>>> pipleline gapA clip 2 <<<<<<<<<<<<< */

    if (sol.config.do_gapA_clip == 1)
      {
       pipelines["pl2/f1/type"] = "clip";
       conduit::Node &clip_params2 = pipelines["pl2/f1/params"];

       clip_params2["topology"] = "topo";
       clip_params2["plane/point/x"] = sol.config.gapA_clip2_point_x;
       clip_params2["plane/point/y"] = sol.config.gapA_clip2_point_y;
       clip_params2["plane/point/z"] = sol.config.gapA_clip2_point_z;
       clip_params2["plane/normal/x"] = sol.config.gapA_clip2_norm_x;
       clip_params2["plane/normal/y"] = sol.config.gapA_clip2_norm_y;
       clip_params2["plane/normal/z"] = sol.config.gapA_clip2_norm_z;

       scenes["s3/plots/p1/type"] = "pseudocolor";
       scenes["s3/plots/p1/pipeline"] = "pl2";
       scenes["s3/plots/p1/field"] = "gapA";

       scenes["s3/plots/p1/min_value"]
	 = matep.gap_A_td(sol.config.Inip, (sol.config.Ttdb0 * matep.Tcp_mK(sol.config.Inip))) * (1. + sol.config.clamp_bias_gapMin);

       scenes["s3/plots/p1/max_value"]
	 = matep.gap_B_td(sol.config.Inip, (sol.config.Ttdb0 * matep.Tcp_mK(sol.config.Inip))) * (1. + sol.config.clamp_bias_gapMax);

       scenes["s3/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s3/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);        
       scenes["s3/renders/r1/image_prefix"] = "insitu/gapA-clip2/gapA-clip2_t-%09d";
       scenes["s3/renders/r1/camera/azimuth"] = sol.config.camera2_azi/*35.0*/;
       scenes["s3/renders/r1/camera/elevation"] = sol.config.camera2_ele/*30.0*/;
      }

    /* >>>>>>>>>>>>>> pipleline gapA slice1 <<<<<<<<<<<<< */
    
    if (sol.config.do_gapA_slice == 1)
      {
       pipelines["pl3/f1/type"] = "exaslice";
       conduit::Node &slice_params1 = pipelines["pl3/f1/params"];

       slice_params1["point/x"] = sol.config.gapA_slice1_point_x;
       slice_params1["point/y"] = sol.config.gapA_slice1_point_y;
       slice_params1["point/z"] = sol.config.gapA_slice1_point_z;
       slice_params1["normal/x"] = sol.config.gapA_slice1_norm_x;
       slice_params1["normal/y"] = sol.config.gapA_slice1_norm_y;
       slice_params1["normal/z"] = sol.config.gapA_slice1_norm_z;

       scenes["s4/plots/p1/type"] = "pseudocolor";
       scenes["s4/plots/p1/pipeline"] = "pl3";
       scenes["s4/plots/p1/field"] = "gapA";

       scenes["s4/plots/p1/min_value"]
	 = 0.0;
    
       scenes["s4/plots/p1/max_value"]
	 = (sol.config.initialConditionT == 2)
	   ? matep.gap_B_td(sol.config.Inip, (sol.config.Ttdb0 * matep.Tcp_mK(sol.config.Inip))) * (1. + sol.config.clamp_bias_gapMax)
	   : matep.gap_B_td(sol.config.Inip, sol.config.IniT) * (1. + sol.config.clamp_bias_gapMax);

       scenes["s4/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s4/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);    
       scenes["s4/renders/r1/image_prefix"] = "insitu/gapA-slice1/gapA-slice1_t-%09d";
       // scenes["s4/renders/r1/camera/azimuth"] = -45./*sol.config.camera1_azi*/;
       // scenes["s4/renders/r1/camera/elevation"] = 0./*sol.config.camera1_ele*/;
      }


    /* >>>>>>>>>>>>>> pipleline gapA slice2 <<<<<<<<<<<<< */
    
    // if (sol.config.do_gapA_slice == 1)
    //   {
    //    pipelines2["pl4/f1/type"] = "slice";
    //    conduit::Node &slice_params2 = pipelines2["pl4/f1/params"];

    //    slice_params2["point/x"] = sol.config.gapA_slice2_point_x;
    //    slice_params2["point/y"] = sol.config.gapA_slice2_point_y;
    //    slice_params2["point/z"] = sol.config.gapA_slice2_point_z;
    //    slice_params2["normal/x"] = sol.config.gapA_slice2_norm_x;
    //    slice_params2["normal/y"] = sol.config.gapA_slice2_norm_y;
    //    slice_params2["normal/z"] = sol.config.gapA_slice2_norm_z;

    //    scenes2["s5/plots/p1/type"] = "pseudocolor";
    //    scenes2["s5/plots/p1/pipeline"] = "pl4";
    //    scenes2["s5/plots/p1/field"] = "gapA";
    //    // scenes["s5/image_prefix"] = "gapA-slice2_t-%09d";       

    //    scenes2["s5/plots/p1/min_value"]
    // 	 = 0.0;
    
    //    scenes2["s5/plots/p1/max_value"]
    // 	 = matep.gap_B_td(sol.config.Inip, (sol.config.Ttdb0 * matep.Tcp_mK(sol.config.Inip))) * (1. + sol.config.clamp_bias_gapMax);

    //    scenes2["s5/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
    //    scenes2["s5/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);    
    //    scenes2["s5/renders/r1/image_prefix"] = "gapA-slice2_t-%09d";
    //    // scenes["s5/renders/r1/camera/azimuth"] = 45. /*sol.config.camera1_azi*/;
    //    // scenes["s5/renders/r1/camera/elevation"] = 0./*sol.config.camera1_ele*/;
    //   }
    
    
    /* >>>>>>>>>>> pipleline isosurfece <<<<<<<<<<<<< */

    if (sol.config.do_gapA_isosurface == 1)
      {
       pipelines["pl5/f1/type"] = "contour";

       conduit::Node &contour_params = pipelines["pl5/f1/params"];
       contour_params["field"] = "gapA";

       //gapA_iso_list_size has to be constexpr in order to fix VLAs warning from Clang
       //defalt set its size to be 2, this should be carefully obeyed in configuration file.       
       const unsigned int iso_list_size = sol.config.iso_values_vector.size();

       //double iso_vals[Temperature_iso_list_size];
       //fix iso_vals array length; fix Clang 19 VLAs warning       
       double iso_vals[2];
       for (unsigned int i = 0; i<iso_list_size; ++i) {iso_vals[i] = sol.config.iso_values_vector[i];}

       contour_params["iso_values"].set(iso_vals, iso_list_size);

       scenes["s6/plots/p1/type"] = "pseudocolor";
       scenes["s6/plots/p1/pipeline"] = "pl5";
       scenes["s6/plots/p1/field"] = "gapA";
       scenes["s6/renders/r1/image_prefix"] = "insitu/gapA-iso/gapA-iso_t-%09d";

       double box_bounds[6] = {0.0, sol.config.lx * sol.config.dx, 0.0, sol.config.ly * sol.config.dx, 0.0, sol.config.lz * sol.config.dx};
       scenes["s6/renders/r1/dataset_bounds"].set(box_bounds,6);
       scenes["s6/renders/r1/render_bg"] = "true";

       scenes["s6/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s6/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);             
       scenes["s6/renders/r1/camera/azimuth"] = sol.config.camera1_azi/*35.0*/;
       scenes["s6/renders/r1/camera/elevation"] = sol.config.camera1_ele/*30.0*/;
    }
    
    /* >>>>>>>>>>>>>> pipleline clip camp <<<<<<<<<<<<< */

    if (sol.config.do_fed_clip == 1)
      {
       pipelines["pl6/f1/type"] = "clip";
       conduit::Node &clip_params3 = pipelines["pl6/f1/params"];

       clip_params3["topology"] = "topo";
       clip_params3["plane/point/x"] = sol.config.fed_clip_point_x;
       clip_params3["plane/point/y"] = sol.config.fed_clip_point_y;
       clip_params3["plane/point/z"] = sol.config.fed_clip_point_z;
       clip_params3["plane/normal/x"] = sol.config.fed_clip_norm_x;
       clip_params3["plane/normal/y"] = sol.config.fed_clip_norm_y;
       clip_params3["plane/normal/z"] = sol.config.fed_clip_norm_z;

       scenes["s7/plots/p1/type"] = "pseudocolor";
       scenes["s7/plots/p1/pipeline"] = "pl6";
       scenes["s7/plots/p1/field"] = "feDensity";

       scenes["s7/plots/p1/min_value"] = matep.f_B_td(sol.config.Inip, sol.config.Ttd_Qend * matep.Tcp_mK(sol.config.Inip)) * (1. + sol.config.clamp_bias_fed_Min);

       scenes["s7/plots/p1/max_value"] = 1.0 * (1. + sol.config.clamp_bias_fed_Max);    

       scenes["s7/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s7/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);        
       scenes["s7/renders/r1/image_prefix"] = "insitu/feDensity-clip_Camp/feDensity-clip_Camp-t-%09d";
       scenes["s7/renders/r1/camera/azimuth"] = sol.config.camera1_azi/*35.0*/;
       scenes["s7/renders/r1/camera/elevation"] = sol.config.camera1_ele/*30.0*/;
      }

    /* >>>>>>>>>>>>>> pipleline Temperaure clip no camp <<<<<<<<<<<<< */

    if (sol.config.do_Temperature_clip == 1)
      {
       pipelines["pl7/f1/type"] = "clip";
       conduit::Node &clip_params4 = pipelines["pl7/f1/params"];

       clip_params4["topology"] = "topo";
       clip_params4["plane/point/x"] = sol.config.Temperature_clip_point_x;
       clip_params4["plane/point/y"] = sol.config.Temperature_clip_point_y;
       clip_params4["plane/point/z"] = sol.config.Temperature_clip_point_z;
       clip_params4["plane/normal/x"] = sol.config.Temperature_clip_norm_x;
       clip_params4["plane/normal/y"] = sol.config.Temperature_clip_norm_y;
       clip_params4["plane/normal/z"] = sol.config.Temperature_clip_norm_z;

       scenes["s8/plots/p1/type"] = "pseudocolor";
       scenes["s8/plots/p1/pipeline"] = "pl7";
       scenes["s8/plots/p1/field"] = "Temperature";

       scenes["s8/plots/p1/min_value"] = sol.config.Ttdb0 * matep.Tcp_mK(sol.config.Inip);

       // this upper limit needs to be calculated carefully
       scenes["s8/plots/p1/max_value"] = sol.config.Temperature_clamp * matep.Tcp_mK(sol.config.Inip);

       scenes["s8/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s8/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);        
       scenes["s8/renders/r1/image_prefix"] = "insitu/Temeperature-clip/Temeperature-clip-t-%09d";
       scenes["s8/renders/r1/camera/azimuth"] = sol.config.camera1_azi/*35.0*/;
       scenes["s8/renders/r1/camera/elevation"] = sol.config.camera1_ele/*30.0*/;
      }

    /* >>>>>>>>>>>>>> pipleline Temperature slice <<<<<<<<<<<<< */
    
    if (sol.config.do_Temperature_slice == 1)
      {
       pipelines2["pl4/f1/type"] = "slice";
       conduit::Node &slice_params2 = pipelines2["pl4/f1/params"];

       slice_params2["point/x"] = sol.config.Temperature_slice_point_x;
       slice_params2["point/y"] = sol.config.Temperature_slice_point_y;
       slice_params2["point/z"] = sol.config.Temperature_slice_point_z;
       slice_params2["normal/x"] = sol.config.Temperature_slice_norm_x;
       slice_params2["normal/y"] = sol.config.Temperature_slice_norm_y;
       slice_params2["normal/z"] = sol.config.Temperature_slice_norm_z;

       scenes2["s5/plots/p1/type"] = "pseudocolor";
       scenes2["s5/plots/p1/pipeline"] = "pl4";
       scenes2["s5/plots/p1/field"] = "Temperature";

       scenes2["s5/plots/p1/min_value"]
	 = sol.config.Ttdb0 * matep.Tcp_mK(sol.config.Inip);
    
       scenes2["s5/plots/p1/max_value"]
	 = matep.Tcp_mK(sol.config.Inip);

       scenes2["s5/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes2["s5/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);    
       scenes2["s5/renders/r1/image_prefix"] = "insitu/Temperature-slice/Temperature-slice_t-%09d";
      }

    
    /* >>>>>>>>>>> Tmeperature isosurfece <<<<<<<<<<<<< */

    if (sol.config.do_Temperature_isosurface == 1)
      {
       pipelines["pl8/f1/type"] = "contour";

       conduit::Node &contour_params = pipelines["pl8/f1/params"];
       contour_params["field"] = "Temperature";

       //Temperature_iso_list_size has to be constexpr in order to fix VLAs warning from Clang
       //defalt set its size to be 2, this should be carefully obeyed in configuration file.
       const unsigned int Temperature_iso_list_size = sol.config.Temperature_iso_values_vector.size();
       //double iso_vals[Temperature_iso_list_size];
       //fix iso_vals array length; fix Clang 19 VLAs warning
       double iso_vals[2];
       for (unsigned int i = 0; i<Temperature_iso_list_size; ++i) {iso_vals[i] = sol.config.Temperature_iso_values_vector[i];}

       contour_params["iso_values"].set(iso_vals, Temperature_iso_list_size);

       scenes["s9/plots/p1/type"] = "pseudocolor";
       scenes["s9/plots/p1/pipeline"] = "pl8";
       scenes["s9/plots/p1/field"] = "Temperature";
       scenes["s9/renders/r1/image_prefix"] = "insitu/Temperature-iso/Temperature-iso_t-%09d";

       double box_bounds[6] = {0.0, sol.config.lx * sol.config.dx, 0.0, sol.config.ly * sol.config.dx, 0.0, sol.config.lz * sol.config.dx};
       scenes["s9/renders/r1/dataset_bounds"].set(box_bounds,6);
       scenes["s9/renders/r1/render_bg"] = "true";

       scenes["s9/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s9/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);             
       scenes["s9/renders/r1/camera/azimuth"] = sol.config.camera1_azi/*35.0*/;
       scenes["s9/renders/r1/camera/elevation"] = sol.config.camera1_ele/*30.0*/;
    }
    
    /* >>>>>>>>>>>>>> pipleline phaseMarker slice <<<<<<<<<<<<< */
    
    if (sol.config.do_phaseMarker_slice == 1)
      {
       pipelines2["pl9/f1/type"] = "exaslice";
       conduit::Node &slice_params3 = pipelines2["pl9/f1/params"];

       slice_params3["point/x"] = sol.config.pMarker_slice_point_x;
       slice_params3["point/y"] = sol.config.pMarker_slice_point_y;
       slice_params3["point/z"] = sol.config.pMarker_slice_point_z;
       slice_params3["normal/x"] = sol.config.pMarker_slice_norm_x;
       slice_params3["normal/y"] = sol.config.pMarker_slice_norm_y;
       slice_params3["normal/z"] = sol.config.pMarker_slice_norm_z;

       scenes["s10/plots/p1/type"] = "pseudocolor";
       scenes["s10/plots/p1/pipeline"] = "pl9";
       scenes["s10/plots/p1/field"] = "phaseMarker";
       scenes["s10/plots/p1/color_table/name"] = "Jet";
       // scenes["s10/plots/p1/color_table/discrete"] = "true";

       //???????????????

       scenes["s10/plots/p1/min_value"]
	 = 1.0f;
    
       scenes["s10/plots/p1/max_value"]
	 = 9.0f; 

       scenes["s10/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s10/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);    
       scenes["s10/renders/r1/image_prefix"] = "insitu/pMarker-Slice/pMarker-slice_t-%09d";
      }

    /* >>>>>>>>>>> phaseMarker isosurfece <<<<<<<<<<<<< */
    
    if (sol.config.do_phaseMarker_isosurface == 1)
      {
       pipelines2["pl10/f1/type"] = "contour";

       conduit::Node &contour_params = pipelines2["pl10/f1/params"];
       contour_params["field"] = "phaseMarker";

       //phaseMarker_iso_list_size has to be constexpr in order to fix VLAs warning from Clang
       //defalt set its size to be 2, this should be carefully obeyed in configuration file.       
       const unsigned int phaseMarker_iso_list_size = sol.config.phaseMarker_iso_values_vector.size();
       //double iso_vals[Temperature_iso_list_size];
       //fix iso_vals array length; fix Clang 19 VLAs warning       
       double iso_vals[2];
       for (unsigned int i = 0; i<phaseMarker_iso_list_size; ++i) {iso_vals[i] = sol.config.phaseMarker_iso_values_vector[i];}

       contour_params["iso_values"].set(iso_vals, phaseMarker_iso_list_size);

       scenes2["s11/plots/p1/type"] = "pseudocolor";
       scenes2["s11/plots/p1/pipeline"] = "pl10";
       scenes2["s11/plots/p1/field"] = "phaseMarker";
       scenes2["s11/plots/p1/color_table/name"] = "Jet";       
       scenes2["s11/renders/r1/image_prefix"] = "insitu/pMarker-iso/pMarker-iso_t-%09d";
       scenes2["s11/renders/r1/screen_annotations"] = (sol.config.remove_screen_annotations == 1) ? "false" : "true";       

       double box_bounds[6] = {0.0, sol.config.lx * sol.config.dx, 0.0, sol.config.ly * sol.config.dx, 0.0, sol.config.lz * sol.config.dx};
       scenes2["s11/renders/r1/dataset_bounds"].set(box_bounds,6);
       scenes2["s11/renders/r1/render_bg"] = "true";

       scenes2["s11/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes2["s11/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);             
       scenes2["s11/renders/r1/camera/azimuth"] = sol.config.camera1_azi/*35.0*/;
       scenes2["s11/renders/r1/camera/elevation"] = sol.config.camera1_ele/*30.0*/;
    }

    /* >>>>>>>>>>> phaseMarker pMarker < 2 fieldclip <<<<<<<<<<<<< */
    
    if (sol.config.do_phaseMarker_fieldclip == 1)
      {
       pipelines2["pl11/f1/type"] = "clip_with_field";

       conduit::Node &clip_params = pipelines2["pl11/f1/params"];
       clip_params["field"] = "phaseMarker";
       clip_params["invert"] = "true";
       clip_params["clip_value"] = 2.;

       scenes2["s12/plots/p1/type"] = "pseudocolor";
       scenes2["s12/plots/p1/pipeline"] = "pl11";
       scenes2["s12/plots/p1/field"] = "phaseMarker";
       scenes2["s12/plots/p1/color_table/name"] = "Default";
       //scenes2["s12/plots/p1/color_table/discrete"] = "true";
       scenes2["s12/renders/r1/image_prefix"] = "insitu/pMarker-fieldclip/pMarker-fieldclip_t-%09d";
       scenes2["s12/renders/r1/screen_annotations"] = (sol.config.remove_screen_annotations == 1) ? "false" : "true";       

       double box_bounds[6] = {0.0, sol.config.lx * sol.config.dx, 0.0, sol.config.ly * sol.config.dx, 0.0, sol.config.lz * sol.config.dx};
       scenes2["s12/renders/r1/dataset_bounds"].set(box_bounds,6);
       scenes2["s12/renders/r1/render_bg"] = "true";

       scenes2["s12/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes2["s12/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);             
       scenes2["s12/renders/r1/camera/azimuth"] = sol.config.camera1_azi/*35.0*/;
       scenes2["s12/renders/r1/camera/elevation"] = sol.config.camera1_ele/*30.0*/;
    }

    /* >>>>>>>>>>> phaseMarker B-isoVolume <<<<<<<<<<<<< */
    
    if (sol.config.do_phaseMarker_fieldclip_Bphase == 1)
      {
       pipelines2["pl12/f1/type"] = "isovolume";

       conduit::Node &clip_params = pipelines2["pl12/f1/params"];
       clip_params["field"] = "phaseMarker";
       clip_params["min_value"] = 3.0;
       clip_params["max_value"] = 6.0;

       scenes2["s13/plots/p1/type"] = "pseudocolor";
       scenes2["s13/plots/p1/pipeline"] = "pl12";
       scenes2["s13/plots/p1/field"] = "phaseMarker";
       scenes2["s13/plots/p1/color_table/name"] = "Green";
       //scenes2["s13/plots/p1/color_table/discrete"] = "true";
       scenes2["s13/renders/r1/image_prefix"] = "insitu/pMarker-isoVolume-Bphase/pMarker-isoVolume-Bphase_t-%09d";
       scenes2["s13/renders/r1/screen_annotations"] = (sol.config.remove_screen_annotations == 1) ? "false" : "true";       

       double box_bounds[6] = {0.0, sol.config.lx * sol.config.dx, 0.0, sol.config.ly * sol.config.dx, 0.0, sol.config.lz * sol.config.dx};
       scenes2["s13/renders/r1/dataset_bounds"].set(box_bounds,6);
       scenes2["s13/renders/r1/render_bg"] = "true";

       scenes2["s13/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes2["s13/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);             
       scenes2["s13/renders/r1/camera/azimuth"] = sol.config.camera1_azi/*35.0*/;
       scenes2["s13/renders/r1/camera/elevation"] = sol.config.camera1_ele/*30.0*/;
    }

    /* >>>>>>>>>>> phaseMarker A-phase fieldclip <<<<<<<<<<<<< */
    
    if (sol.config.do_phaseMarker_fieldclip_Aphase == 1)
      {
       pipelines2["pl13/f1/type"] = "clip_with_field";

       conduit::Node &clip_params = pipelines2["pl13/f1/params"];
       clip_params["field"] = "phaseMarker";
       clip_params["invert"] = "false";
       clip_params["clip_value"] = 8.;

       scenes2["s14/plots/p1/type"] = "pseudocolor";
       scenes2["s14/plots/p1/pipeline"] = "pl13";
       scenes2["s14/plots/p1/field"] = "phaseMarker";
       scenes2["s14/plots/p1/color_table/name"] = "Cold and Hot";
       //scenes2["s14/plots/p1/color_table/discrete"] = "true";
       scenes2["s14/renders/r1/image_prefix"] = "insitu/pMarker-fieldclip-Aphase/pMarker-fieldclip-Aphase_t-%09d";
       scenes2["s14/renders/r1/screen_annotations"] = (sol.config.remove_screen_annotations == 1) ? "false" : "true";

       double box_bounds[6] = {0.0, sol.config.lx * sol.config.dx, 0.0, sol.config.ly * sol.config.dx, 0.0, sol.config.lz * sol.config.dx};
       scenes2["s14/renders/r1/dataset_bounds"].set(box_bounds,6);
       scenes2["s14/renders/r1/render_bg"] = "true";

       scenes2["s14/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes2["s14/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);             
       scenes2["s14/renders/r1/camera/azimuth"] = sol.config.camera1_azi/*35.0*/;
       scenes2["s14/renders/r1/camera/elevation"] = sol.config.camera1_ele/*30.0*/;
      }
    
    /* >>>>>>>>>>>>>> pipleline U(1) 3 \phi slice <<<<<<<<<<<<< */
    
    if (sol.config.do_U13phi_slice == 1)
      {
       pipelines2["pl14/f1/type"] = "exaslice";
       conduit::Node &slice_params4 = pipelines2["pl14/f1/params"];

       slice_params4["point/x"] = sol.config.U13phi_slice_point_x;
       slice_params4["point/y"] = sol.config.U13phi_slice_point_y;
       slice_params4["point/z"] = sol.config.U13phi_slice_point_z;
       slice_params4["normal/x"] = sol.config.U13phi_slice_norm_x;
       slice_params4["normal/y"] = sol.config.U13phi_slice_norm_y;
       slice_params4["normal/z"] = sol.config.U13phi_slice_norm_z;

       /* ----------------------------------------------------------
        * custom hsv style Cyclic color map, 
        * RGBA control points are defined in hsvCyclicSchemeRGBA.hpp
        * ----------------------------------------------------------
        */
       conduit::Node hsv_Cyclic;
       hsv_Cyclic["r"].set_external(hsv_r);
       hsv_Cyclic["g"].set_external(hsv_g);
       hsv_Cyclic["b"].set_external(hsv_b);
       hsv_Cyclic["a"].set_external(hsv_a);
       hsv_Cyclic["position"].set_external(hsv_value_point_position);
       
       scenes["s15/plots/p1/type"] = "pseudocolor";
       scenes["s15/plots/p1/pipeline"] = "pl14";
       scenes["s15/plots/p1/field"] = "U1_3phi";
       scenes["s15/plots/p1/min_value"] = -3.141592653589793;
       scenes["s15/plots/p1/max_value"] = 3.141592653589793;
       //scenes["s15/plots/p1/color_table/name"] = "Blue to Orange";
       scenes["s15/plots/p1/color_table/control_points"] = hsv_Cyclic;
       // scenes["s15/plots/p1/color_table/discrete"] = "true";

       scenes["s15/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s15/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);    
       scenes["s15/renders/r1/image_prefix"] = "insitu/U1_3phi-Slice/U1_3phi-slice_t-%09d";
      }

    /* >>>>>>>>>>>>>> pipleline l^2 slice <<<<<<<<<<<<< */
    
    if (sol.config.do_l_sq_slice == 1)
      {
       pipelines2["pl15/f1/type"] = "exaslice";
       conduit::Node &slice_params5 = pipelines2["pl15/f1/params"];

       slice_params5["point/x"] = sol.config.l_sq_slice_point_x;
       slice_params5["point/y"] = sol.config.l_sq_slice_point_y;
       slice_params5["point/z"] = sol.config.l_sq_slice_point_z;
       slice_params5["normal/x"] = sol.config.l_sq_slice_norm_x;
       slice_params5["normal/y"] = sol.config.l_sq_slice_norm_y;
       slice_params5["normal/z"] = sol.config.l_sq_slice_norm_z;

       scenes["s16/plots/p1/type"] = "pseudocolor";
       scenes["s16/plots/p1/pipeline"] = "pl15";
       scenes["s16/plots/p1/field"] = "l_Sq";
       scenes["s16/plots/p1/color_table/name"] = "Cool to Warm Extended";
       // scenes["s16/plots/p1/color_table/discrete"] = "true";

       scenes["s16/plots/p1/min_value"]
	 = 0.f; /*1.f - sol.config.lVec_SqlTol;*/
    
       scenes["s16/plots/p1/max_value"]
	 = 1.f + sol.config.lVec_SqlTol; 

       scenes["s16/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s16/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);    
       scenes["s16/renders/r1/image_prefix"] = "insitu/lVec_norm2-Slice/lVec_norm2-slice_t-%09d";
      }

    
} // defineActions() call end here


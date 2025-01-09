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
       scenes["s2/renders/r1/image_prefix"] = "gapA-clip1_t-%05d";
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
       scenes["s3/renders/r1/image_prefix"] = "gapA-clip2-t-%05d";
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
       // scenes["s4/image_prefix"] = "gapA-slice1_t-%05d";       

       scenes["s4/plots/p1/min_value"]
	 = 0.0;
    
       scenes["s4/plots/p1/max_value"]
	 = matep.gap_B_td(sol.config.Inip, (sol.config.Ttdb0 * matep.Tcp_mK(sol.config.Inip))) * (1. + sol.config.clamp_bias_gapMax);

       scenes["s4/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s4/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);    
       scenes["s4/renders/r1/image_prefix"] = "gapA-slice1_t-%05d";
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
    //    // scenes["s5/image_prefix"] = "gapA-slice2_t-%05d";       

    //    scenes2["s5/plots/p1/min_value"]
    // 	 = 0.0;
    
    //    scenes2["s5/plots/p1/max_value"]
    // 	 = matep.gap_B_td(sol.config.Inip, (sol.config.Ttdb0 * matep.Tcp_mK(sol.config.Inip))) * (1. + sol.config.clamp_bias_gapMax);

    //    scenes2["s5/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
    //    scenes2["s5/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);    
    //    scenes2["s5/renders/r1/image_prefix"] = "gapA-slice2_t-%05d";
    //    // scenes["s5/renders/r1/camera/azimuth"] = 45. /*sol.config.camera1_azi*/;
    //    // scenes["s5/renders/r1/camera/elevation"] = 0./*sol.config.camera1_ele*/;
    //   }
    
    
    /* >>>>>>>>>>> pipleline isosurfece <<<<<<<<<<<<< */

    if (sol.config.do_gapA_isosurface == 1)
      {
       pipelines["pl5/f1/type"] = "contour";

       conduit::Node &contour_params = pipelines["pl5/f1/params"];
       contour_params["field"] = "gapA";

       const unsigned int iso_list_size = sol.config.iso_values_vector.size();
       double iso_vals[iso_list_size];
       for (unsigned int i = 0; i<iso_list_size; ++i) {iso_vals[i] = sol.config.iso_values_vector[i];}

       contour_params["iso_values"].set(iso_vals, iso_list_size);

       scenes["s6/plots/p1/type"] = "pseudocolor";
       scenes["s6/plots/p1/pipeline"] = "pl5";
       scenes["s6/plots/p1/field"] = "gapA";
       scenes["s6/renders/r1/image_prefix"] = "gapA-iso_t-%05d";

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
       scenes["s7/renders/r1/image_prefix"] = "feDensity-clip_Camp-t-%05d";
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
       scenes["s8/renders/r1/image_prefix"] = "Temeperature-clip-t-%05d";
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
       scenes2["s5/renders/r1/image_prefix"] = "Temperature-slice_t-%05d";
      }

    
    /* >>>>>>>>>>> Tmeperature isosurfece <<<<<<<<<<<<< */

    if (sol.config.do_Temperature_isosurface == 1)
      {
       pipelines["pl8/f1/type"] = "contour";

       conduit::Node &contour_params = pipelines["pl8/f1/params"];
       contour_params["field"] = "Temperature";

       const unsigned int Temperature_iso_list_size = sol.config.Temperature_iso_values_vector.size();
       double iso_vals[Temperature_iso_list_size];
       for (unsigned int i = 0; i<Temperature_iso_list_size; ++i) {iso_vals[i] = sol.config.Temperature_iso_values_vector[i];}

       contour_params["iso_values"].set(iso_vals, Temperature_iso_list_size);

       scenes["s9/plots/p1/type"] = "pseudocolor";
       scenes["s9/plots/p1/pipeline"] = "pl8";
       scenes["s9/plots/p1/field"] = "Temperature";
       scenes["s9/renders/r1/image_prefix"] = "Temperature-iso_t-%05d";

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
       pipelines2["pl5/f1/type"] = "exaslice";
       conduit::Node &slice_params3 = pipelines2["pl5/f1/params"];

       slice_params3["point/x"] = sol.config.pMarker_slice_point_x;
       slice_params3["point/y"] = sol.config.pMarker_slice_point_y;
       slice_params3["point/z"] = sol.config.pMarker_slice_point_z;
       slice_params3["normal/x"] = sol.config.pMarker_slice_norm_x;
       slice_params3["normal/y"] = sol.config.pMarker_slice_norm_y;
       slice_params3["normal/z"] = sol.config.pMarker_slice_norm_z;

       scenes["s6/plots/p1/type"] = "pseudocolor";
       scenes["s6/plots/p1/pipeline"] = "pl5";
       scenes["s6/plots/p1/field"] = "phaseMarker";
       scenes["s6/plots/p1/color_table/name"] = "Jet";

       //???????????????

       scenes["s6/plots/p1/min_value"]
	 = 0.0f;/*sol.config.Ttdb0 * matep.Tcp_mK(sol.config.Inip);*/
    
       scenes["s6/plots/p1/max_value"]
	 = 4.0f; /*matep.Tcp_mK(sol.config.Inip);*/

       scenes["s6/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
       scenes["s6/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);    
       scenes["s6/renders/r1/image_prefix"] = "pMarker-slice_t-%05d";
      }
    
} // defineActions() call end here


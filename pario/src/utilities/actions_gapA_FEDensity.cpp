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


void parIO::defineActions_gapA_FEDensity(glsol &sol) {

    conduit::Node &add_act = actions.append();    
    add_act["action"] = "add_scenes";
    conduit::Node &scenes = add_act["scenes"];

    // background color vecctor
    double bg_colvec[3] = {0.8, 0.8, 0.8};
    // foreground color vector
    double fg_colvec[3] = {0., 0., 0.};
    
    /* >>>>>>>>>>>>>>   gapA 1st scene    <<<<<<<<<<<<< */
    scenes["s1/plots/p1/type"] = "pseudocolor";
    scenes["s1/plots/p1/field"] = "gapAOrdered";

    // color map clamping. min_value will be set to 0.0 if initialCondtion is 2, 7 i.e., normal_phase, Aphase_full
    scenes["s1/plots/p1/min_value"]
      = (sol.config.initialCondition == 7 || sol.config.initialCondition == 2 || sol.config.initialCondition == 0)
         ? 0.0 :  (matep.gap_A_td(sol.config.Inip, sol.config.IniT)
		   * (1. + sol.config.clamp_bias_gapMin));

    scenes["s1/plots/p1/max_value"]
      = matep.gap_B_td(sol.config.Inip, sol.config.IniT)
        * (1. + sol.config.clamp_bias_gapMax);
    
    scenes["s1/renders/r1/image_prefix"] = "gapA_t-%04d";
    scenes["s1/renders/r1/camera/azimuth"] = sol.config.camera_azi/*35.0*/;
    scenes["s1/renders/r1/camera/elevation"] = sol.config.camera_ele/*30.0*/;

    /* >>>>>>>>>>>>>> pipleline clip <<<<<<<<<<<<< */

    conduit::Node &add_act2 = actions.append();
    add_act2["action"] = "add_pipelines";
    conduit::Node &pipelines = add_act2["pipelines"];

    pipelines["pl1/f1/type"] = "clip";
    conduit::Node &clip_params = pipelines["pl1/f1/params"];

    clip_params["topology"] = "topo";
    clip_params["plane/point/x"] = sol.config.gapA_clip_point_x;
    clip_params["plane/point/y"] = sol.config.gapA_clip_point_y;
    clip_params["plane/point/z"] = sol.config.gapA_clip_point_z;
    clip_params["plane/normal/x"] = sol.config.gapA_clip_norm_x;
    clip_params["plane/normal/y"] = sol.config.gapA_clip_norm_y;
    clip_params["plane/normal/z"] = sol.config.gapA_clip_norm_z;

    scenes["s2/plots/p1/type"] = "pseudocolor";
    scenes["s2/plots/p1/pipeline"] = "pl1";
    scenes["s2/plots/p1/field"] = "gapAOrdered";

    scenes["s2/plots/p1/min_value"] = (sol.config.initialCondition == 2 || sol.config.initialCondition == 0)
                                       ? 0.0 :  matep.gap_A_td(sol.config.Inip, sol.config.IniT) * (1. + sol.config.clamp_bias_gapMin);
    scenes["s2/plots/p1/max_value"] = matep.gap_B_td(sol.config.Inip, (sol.config.Ttd_Qend * matep.Tcp_mK(sol.config.Inip))) * (1 + sol.config.clamp_bias_gapMax);

    scenes["s2/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
    scenes["s2/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);    
    scenes["s2/renders/r1/image_prefix"] = "gapA-clip_t-%04d";
    scenes["s2/renders/r1/camera/azimuth"] = sol.config.camera_azi/*35.0*/;
    scenes["s2/renders/r1/camera/elevation"] = sol.config.camera_ele/*30.0*/;

    /* >>>>>>>>>>>>>> pipleline clip no camp <<<<<<<<<<<<< */

    // conduit::Node &add_act3 = actions.append();
    // add_act3["action"] = "add_pipelines";
    // conduit::Node &pipelines2 = add_act3["pipelines"];

    pipelines["pl2/f1/type"] = "clip";
    conduit::Node &clip_params2 = pipelines["pl2/f1/params"];

    clip_params2["topology"] = "topo";
    clip_params2["plane/point/x"] = sol.config.gapA_clip_point_x;
    clip_params2["plane/point/y"] = sol.config.gapA_clip_point_y;
    clip_params2["plane/point/z"] = sol.config.gapA_clip_point_z;
    clip_params2["plane/normal/x"] = sol.config.gapA_clip_norm_x;
    clip_params2["plane/normal/y"] = sol.config.gapA_clip_norm_y;
    clip_params2["plane/normal/z"] = sol.config.gapA_clip_norm_z;

    scenes["s3/plots/p1/type"] = "pseudocolor";
    scenes["s3/plots/p1/pipeline"] = "pl2";
    scenes["s3/plots/p1/field"] = "gapAOrdered";

    scenes["s3/plots/p1/min_value"] = 0.0;
    // scenes["s3/plots/p1/max_value"] = matep.gap_B_td(sol.config.Inip, sol.config.IniT) * (1. + sol.config.clamp_bias_gapMax);
    scenes["s3/plots/p1/max_value"] = matep.gap_B_td(sol.config.Inip, (sol.config.Ttd_Qend * matep.Tcp_mK(sol.config.Inip))) * (1. + sol.config.clamp_bias_gapMax);    

    scenes["s3/renders/r1/bg_color"].set_float64_ptr(bg_colvec, 3);
    scenes["s3/renders/r1/fg_color"].set_float64_ptr(fg_colvec, 3);        
    scenes["s3/renders/r1/image_prefix"] = "gapA-clip_noCamp-t-%04d";
    scenes["s3/renders/r1/camera/azimuth"] = sol.config.camera_azi/*35.0*/;
    scenes["s3/renders/r1/camera/elevation"] = sol.config.camera_ele/*30.0*/;
    
          
    /* >>>>>>>>>>>>>>  fe_density  2nd scene    <<<<<<<<<<<<< */    
    scenes["s5/plots/p1/type"] = "pseudocolor";
    scenes["s5/plots/p1/field"] = "feDensityOrdered";

    scenes["s5/plots/p1/min_value"] = sol.config.clamp_fed_Min;					  
    scenes["s5/plots/p1/max_value"] = sol.config.clamp_fed_Max;
    
    scenes["s5/renders/r1/image_prefix"] = "feDensity_t-%04d";
    scenes["s5/renders/r1/camera/azimuth"] = sol.config.camera_azi/*35.0*/;
    scenes["s5/renders/r1/camera/elevation"] = sol.config.camera_ele/*30.0*/;
       
} // defineActions() call end here


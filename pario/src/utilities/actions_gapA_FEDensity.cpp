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

    /* >>>>>>>>>>>>>>   1st scene    <<<<<<<<<<<<< */
    scenes["s1/plots/p1/type"] = "pseudocolor";
    scenes["s1/plots/p1/field"] = "gapAOrdered";

    // color map clamping. min_value will be set to 0.0 if initialCondtion is 2 i.e., normal_phase
    scenes["s1/plots/p1/min_value"]
      = (sol.config.initialCondition == 2 || sol.config.initialCondition == 0)
         ? 0.0 :  (matep.gap_A_td(sol.config.Inip, sol.config.IniT)
		   * (1. + sol.config.clamp_bias_gapMin));

    scenes["s1/plots/p1/max_value"]
      = matep.gap_B_td(sol.config.Inip, sol.config.IniT)
        * (1. + sol.config.clamp_bias_gapMax);
    
    scenes["s1/renders/r1/image_prefix"] = "gapA_t-%04d";
    scenes["s1/renders/r1/camera/azimuth"] = sol.config.camera_azi/*35.0*/;
    scenes["s1/renders/r1/camera/elevation"] = sol.config.camera_ele/*30.0*/;
      
    // /* >>>>>>>>>>>>>>   2nd scene    <<<<<<<<<<<<< */    
    scenes["s5/plots/p1/type"] = "pseudocolor";
    scenes["s5/plots/p1/field"] = "feDensityOrdered";

    scenes["s5/plots/p1/min_value"] = 0.0;					  
    scenes["s5/plots/p1/max_value"] = 5.0;
    
    scenes["s5/renders/r1/image_prefix"] = "feDensity_t-%04d";
    scenes["s5/renders/r1/camera/azimuth"] = sol.config.camera_azi/*35.0*/;
    scenes["s5/renders/r1/camera/elevation"] = sol.config.camera_ele/*30.0*/;
       
} // defineActions() call end here


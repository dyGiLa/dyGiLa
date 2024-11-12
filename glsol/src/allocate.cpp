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


const std::string glsol::allocate(const std::string &fname, int argc, char **argv) {
  
    hila::initialize(argc, argv);
    hila::input parameters(fname);

    // First set the simulation box
    
    config.lx = parameters.get("Nx");
    config.ly = parameters.get("Ny");
    config.lz = parameters.get("Nz");
    config.dx = parameters.get("dx");
    config.dtdxRatio = parameters.get("dtdxRatio");

    //config.tThermalizationWaiting = parameters.get("tThermalizationWaiting");
    //config.tauQ1                  = parameters.get("tauQ1");
    //config.tauQ2                  = parameters.get("tauQ2");    
    //config.has1stQStop            = parameters.get_item("has1stQStop",{"no", "yes"});
    //if (config.has1stQStop == 1)
    //  {
    //	config.Ttd_Q1st = parameters.get("Ttd_Q1st");
    //	config.tQ1Waiting = parameters.get("tQ1Waiting");
    //  }    
    //config.Ttd_Qend = parameters.get("Ttd_Qend");
    config.tStart = parameters.get("tStart");
    config.tEnd = parameters.get("tEnd");
    config.tdif = parameters.get("tdif");
    config.difFac = parameters.get("difFac");
    config.tdis = parameters.get("tdis");

    /******************************************/
    /*----   gamma as a complex number   -----*/
    /*   and gamma changes at certain stages  */
    /*   gamma1 and gamma2 are used to hold   */
    /*   different values of gamma.           */
    /*   gamma change is triggered in main()  */
    /******************************************/    
    //config.gamma = parameters.get("gamma");
    std::vector<real_t> tmp1 = parameters.get("gamma1");
    std::vector<real_t> tmp2 = parameters.get("gamma2");    

    if (tmp1.size() < 1 || tmp1.size() >2)
      {
	hila::out0 << "error: gamma1 must be initialized by at least one real or imagnary parts" << std::endl;
	hila::error("\n");
      }
    else if (tmp1.size() == 1)
      {
	hila::out0 << " tmp1.size() = " << tmp1.size() << "\n";
        config.gamma1.real() = tmp1[0];
        config.gamma1.imag() = 0.;
      }
    else  // tmp.size() == 2
      {
        hila::out0 << " tmp1.size() = " << tmp1.size() << "\n";
        config.gamma1.real() = tmp1[0];
        config.gamma1.imag() = tmp1[1];
      }

    if (tmp2.size() < 1 || tmp2.size() >2)
      {
	hila::out0 << "error: gamma must be initialized by at least one real or imagnary parts" << std::endl;
	hila::error("\n");
      }
    else if (tmp2.size() == 1)
      {
	hila::out0 << " tmp2.size() = " << tmp2.size() << "\n";
        config.gamma2.real() = tmp2[0];
        config.gamma2.imag() = 0.;
      }
    else  // tmp.size() == 2
      {
        hila::out0 << " tmp2.size() = " << tmp2.size() << "\n";
        config.gamma2.real() = tmp2[0];
        config.gamma2.imag() = tmp2[1];
      }

    // conunter for triggering gamma value change 
    config.gammaoffc = parameters.get("gammaoffc");
    
    /******************************************/
    /*        gamma as complex number         */
    /*  with gamma1(2) reading end here       */
    /******************************************/        

    // Initialize the order parameter
    
    config.initialCondition = parameters.get_item("initialCondition",{"gaussrand"             //0
								      ,"kgaussrand"           //1
								      ,"normal_phase_real1"   //2
								      ,"normal_phase_real2"   //3
								      ,"normal_phase_complex" //4
								      ,"Bphase"               //5
								      ,"Aphase_partial1"      //6
								      ,"Aphase_full"          //7
                                                                      ,"BinA"});              //8

    hila::out0 << "config.initialCondition is "
	       << config.initialCondition
	       << "\n";
    
    config.variance_sigma = parameters.get("sigma");
    
    config.seed = parameters.get("seed");
    config.IniMod = parameters.get("IniMod");
    config.Inilc = parameters.get("Inilc");

    //Boundary conditions for order parameter

    config.bcx0 = parameters.get_item("bcx0",{"periodic", "A", "B","Normal", "Maximal", "Robin"});
    config.bcxN = parameters.get_item("bcxN",{"periodic", "A", "B","Normal", "Maximal", "Robin"});
    config.bcy0 = parameters.get_item("bcy0",{"periodic", "A", "B","Normal", "Maximal", "Robin"});
    config.bcyN = parameters.get_item("bcyN",{"periodic", "A", "B","Normal", "Maximal", "Robin"});
    config.bcz0 = parameters.get_item("bcz0",{"periodic", "A", "B","Normal", "Maximal", "Robin"});
    config.bczN = parameters.get_item("bczN",{"periodic", "A", "B","Normal", "Maximal", "Robin"});
    config.robinbt_x0 = parameters.get("robinbt_x0");
    config.robinbt_xN = parameters.get("robinbt_xN");
    config.robinbt_y0 = parameters.get("robinbt_y0");
    config.robinbt_yN = parameters.get("robinbt_yN");
    config.robinbt_z0 = parameters.get("robinbt_z0");
    config.robinbt_zN = parameters.get("robinbt_zN");

    //set how the order parameter is going to be evolved
    
    config.useTbath = parameters.get_item("useTbath",{"no","yes"});
    config.Tbath_start = parameters.get("Tbath_start");

    //set boundary conditions for temperature
    
    // Initialize temperature
    
    config.initialConditionT = parameters.get_item("initialConditionT",{"constant","sine","hotspot"});
    if(config.initialConditionT == 0)
      {
	config.IniT = parameters.get("IniT");
      }
    else if (config.initialConditionT == 1)
      {
	config.IniT = parameters.get("IniT");
	config.ampT = parameters.get("ampT");
      }
    else if (config.initialConditionT == 2)
      {
	config.IniT = parameters.get("IniT");
        config.ampT = parameters.get("ampT");
	config.sigTx = parameters.get("sigTx");
	config.sigTy = parameters.get("sigTy");
	config.sigTz = parameters.get("sigTz");
      }

    // set how the temperature is going to be evolved

    config.evolveT = parameters.get_item("evolveT",{"no","yes"});
#ifndef T_FIELD
    hila::out0 << "Tenperature set to a scalar, same over the whole box, only manual T evolution available"<< "\n";
    //config.evolveT=0;                                                                                                                                                                                                                                                                    
#endif

    if(config.evolveT ==1)
      {
        config.Tevolvetype = parameters.get_item("Tevolvetype",{"heat","wave","homogeneousQuench"});
	config.start_evolveT = parameters.get("start_evolveT");
        if (config.Tevolvetype == 0 || config.Tevolvetype == 1)
          {
           config.startdiffT = parameters.get("startdiffT");
           config.diffT = parameters.get("diffT");
          }
      }

    config.bcTx0 = parameters.get_item("bcTx0",{"periodic", "fixed"});
    config.bcTxN = parameters.get_item("bcTxN",{"periodic", "fixed"});
    config.bcTy0 = parameters.get_item("bcTy0",{"periodic", "fixed"});
    config.bcTyN = parameters.get_item("bcTyN",{"periodic", "fixed"});
    config.bcTz0 = parameters.get_item("bcTz0",{"periodic", "fixed"});
    config.bcTzN = parameters.get_item("bcTzN",{"periodic", "fixed"});

    
    // Initialize pressure
    config.initialConditionp = parameters.get_item("initialConditionp",{"constant"});
    config.Inip	= parameters.get("Inip");
    

    //Initialize magnetic field
    config.IniHx = parameters.get("IniHx");
    config.IniHy = parameters.get("IniHy");
    config.IniHz = parameters.get("IniHz");

    config.tStats = parameters.get("tStats");
    config.nOutputs = parameters.get("nOutputs");

    // output_file is the saving path of output file, which offered in congigration file
    const std::string output_file = parameters.get("output_file");

    // xdmf file name, which is provided through config file
    config.xmf2_fname = parameters.get("xmf2_file");


        
    /*----------------------------------------*/
    /* Parallel IO Engine control parameters  */
    /*----------------------------------------*/
    config.hdf5_A_matrix_output        = parameters.get_item("hdf5_A_matrix_output",{"no","yes"});
    if (config.hdf5_A_matrix_output == 1)
      {
       config.hdf5Ststart = parameters.get("hdf5Ststart");
       config.hdf5Stend   = parameters.get("hdf5Stend");       
      }
    
    // config.hdf5_trA_output             = parameters.get_item("hdf5_trA_output",{"no","yes"});
    // config.hdf5_eigvA_output           = parameters.get_item("hdf5_eigvA_output",{"no","yes"});
    config.hdf5_mass_current_output    = parameters.get_item("hdf5_mass_current_output",{"no","yes"});
    config.hdf5_spin_current_output    = parameters.get_item("hdf5_spin_current_output",{"no","yes"});        

    config.do_gapA_clip         = parameters.get_item("do_gapA_clip",{"no","yes"});
    if ( config.do_gapA_clip ==1 )
      {
        config.gapA_clip_point_x = parameters.get("gapA_clip_point_x");
	config.gapA_clip_point_y = parameters.get("gapA_clip_point_y");
	config.gapA_clip_point_z = parameters.get("gapA_clip_point_z");
        config.gapA_clip_norm_x = parameters.get("gapA_clip_norm_x");
	config.gapA_clip_norm_y = parameters.get("gapA_clip_norm_y");
	config.gapA_clip_norm_z = parameters.get("gapA_clip_norm_z");	
      } // gapA clip control parammeters

    config.do_fed_clip         = parameters.get_item("do_fed_clip",{"no","yes"});
    if ( config.do_fed_clip ==1 )
      {
        config.fed_clip_point_x = parameters.get("fed_clip_point_x");
	config.fed_clip_point_y = parameters.get("fed_clip_point_y");
	config.fed_clip_point_z = parameters.get("fed_clip_point_z");
        config.fed_clip_norm_x = parameters.get("fed_clip_norm_x");
	config.fed_clip_norm_y = parameters.get("fed_clip_norm_y");
	config.fed_clip_norm_z = parameters.get("fed_clip_norm_z");	
      } // gapA clip control parammeters
    
    
    config.do_gapA_isosurface   = parameters.get_item("do_gapA_isosurface",{"no","yes"});
    if ( config.do_gapA_isosurface == 1 )
      {
        config.iso_values_vector = parameters.get("iso_values_vector");
      } // gapA clip control parammeters

    
    // config.do_gapA_3slice       = parameters.get_item("do_gapA_3slice",{"no","yes"});
    // config.do_fe_slice          = parameters.get_item("do_fe_slice",{"no","yes"});
    // config.do_gapA_slice        = parameters.get_item("do_gapA_slice",{"no","yes"});            
    
    config.clamp_bias_gapMin = parameters.get("clamp_bias_gapMin");
    config.clamp_bias_gapMax = parameters.get("clamp_bias_gapMax");
    config.clamp_bias_fed_Min = parameters.get("clamp_bias_fed_Min");
    config.clamp_bias_fed_Max = parameters.get("clamp_bias_fed_Max");
    
    config.camera_azi = parameters.get("camera_azi");
    config.camera_ele = parameters.get("camera_ele");
    /*----------------------------------------*/
    /* Parallel IO Engine parameters end      */
    /*----------------------------------------*/
	
    config.dt = config.dx * config.dtdxRatio;
    t = config.tStart;

    // setup the hila lattice geometry 
    CoordinateVector box_dimensions = {config.lx, config.ly, config.lz};
    lattice.setup(box_dimensions);
    hila::seed_random(config.seed);

    return output_file;
    
} // allocate() function ends here


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
//#include "matep.hpp"


const std::string glsol::allocate(const std::string &fname, int argc, char **argv) {
  
    hila::initialize(argc, argv);
    hila::input parameters(fname);
    
    config.lx = parameters.get("Nx");
    config.ly = parameters.get("Ny");
    config.lz = parameters.get("Nz");
    config.dx = parameters.get("dx");
    config.dtdxRatio = parameters.get("dtdxRatio");

    /*********************************************/
    /* >>>>>>  Homogenous quench parameters <<<< */    
    /*********************************************/    
    config.tThermalizationWaiting = parameters.get("tThermalizationWaiting");
    config.tauQ1                  = parameters.get("tauQ1");
    config.tauQ2                  = parameters.get("tauQ2");    
    config.has1stQStop            = parameters.get_item("has1stQStop",{"no", "yes"});
    if (config.has1stQStop == 1)
      {
	config.Ttd_Q1st = parameters.get("Ttd_Q1st");
	config.tQ1Waiting = parameters.get("tQ1Waiting");
      }    
    config.Ttd_Qend = parameters.get("Ttd_Qend");

    /*********************************************/
    /* > Homogenous quench parameters end here<< */    
    /*********************************************/    
    
    
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

    // switch for turning on the T-dependent gamma
    config.TDependnetgamma = parameters.get_item("TDependnetgamma",{"no", "yes"});
    
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

    config.initialCondition = parameters.get_item("initialCondition",{"gaussrand"             //0
								      ,"kgaussrand"           //1
								      ,"normal_phase_real1"   //2
								      ,"normal_phase_real2"   //3
								      ,"normal_phase_complex" //4
								      ,"Bphase"               //5
								      ,"Aphase_partial1"      //6
								      ,"Aphase_full"          //7
	                                                              ,"hotblob"});           //8
                                                                       

    hila::out0 << "config.initialCondition is "
	       << config.initialCondition
	       << "\n";
    
    config.variance_sigma = parameters.get("sigma");
    
    config.seed = parameters.get("seed");
    config.IniMod = parameters.get("IniMod");
    config.Inilc = parameters.get("Inilc");

    //initialCondition-T
    config.initialConditionT = parameters.get_item("initialConditionT",{"constant","sine","hotblob"});
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
       /*********************************************/
       /* >>>> spherical hot blob parameters <<<<<  */    
       /*********************************************/    	
        config.Ttdb1         = parameters.get("Ttdb1");
        config.Ttdb0         = parameters.get("Ttdb0");
        config.t1            = parameters.get("t1");
       /*********************************************/
       /* >>>>>> hot bloob parameters end here  <<< */    
       /*********************************************/    	
      }

    //initialCondition-p
    config.initialConditionp = parameters.get_item("initialConditionp",{"constant"});
    if (config.initialConditionp == 0) { config.Inip = parameters.get("Inip"); }

    //initialCondition-Hfield
    config.withHfield = parameters.get_item("withHfield",{"no", "yes"});
    if (config.withHfield == 1)
      {
       config.initialConditionH = parameters.get_item("initialConditionH",{"constant"});
       std::vector<real_t> temp	= parameters.get("InitH");
       foralldir(al){ config.InitH.e(al) = temp[al]; }
      }
    
    
    config.tStats = parameters.get("tStats");
    config.nOutputs = parameters.get("nOutputs");

    // output_file is the saving path of output file, which offered in congigration file
    const std::string output_file = parameters.get("output_file");

    // xdmf file name, which is provided through config file
    config.xmf2_fname = parameters.get("xmf2_file");

    /*----------------------------------------*/
    /* >>>>>>>  boundary conditions  <<<<<<<<<*/
    /*----------------------------------------*/
    config.BCs1 = parameters.get_item("BCs1",{"periodic",
					       "AB",
					       "PairBreaking",
                                               "PB_y",
                                               "PairB_yz",
                                               "BB",
                                               "phaseVortices"});
    
    config.BCs2 = parameters.get_item("BCs2",{"periodic",
					      "AB",
					      "PairBreaking",
                                              "PB_y",
                                              "PairB_yz",
                                              "BB",
                                              "phaseVortices"});
    
    // config.Wn = parameters.get("BoundaryPhaseWindingNO");
    
    config.BCchangec = parameters.get("BCchangec");
       
    /*if(config.positions==1)
      {
	config.npositionout = parameters.get("npositionout");
	config.write_phases = parameters.get_item("write_phases",{"no","yes"});
	config.write_eigen = parameters.get_item("write_eigen",{"no","yes"});
      }*/
        
    config.evolveT = parameters.get_item("evolveT",{"no","yes"});
    if(config.evolveT ==1
       && config.initialConditionT != 2 )
      {
	config.Tevolvetype = parameters.get_item("Tevolvetype",{"heat","wave","homogeneousQuench"});
	if (config.Tevolvetype == 0 || config.Tevolvetype == 1)
	  {
	   config.startdiffT = parameters.get("startdiffT");
	   config.diffT = parameters.get("diffT");
	  }
      }

    config.useTbath = parameters.get_item("useTbath",{"no","yes"});
    config.Tbath_start = parameters.get("Tbath_start");

    config.ptol = parameters.get("ptol");
    
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
        config.gapA_clip1_point_x = parameters.get("gapA_clip1_point_x");
	config.gapA_clip1_point_y = parameters.get("gapA_clip1_point_y");
	config.gapA_clip1_point_z = parameters.get("gapA_clip1_point_z");
        config.gapA_clip1_norm_x = parameters.get("gapA_clip1_norm_x");
	config.gapA_clip1_norm_y = parameters.get("gapA_clip1_norm_y");
	config.gapA_clip1_norm_z = parameters.get("gapA_clip1_norm_z");

        config.gapA_clip2_point_x = parameters.get("gapA_clip2_point_x");
	config.gapA_clip2_point_y = parameters.get("gapA_clip2_point_y");
	config.gapA_clip2_point_z = parameters.get("gapA_clip2_point_z");
        config.gapA_clip2_norm_x = parameters.get("gapA_clip2_norm_x");
	config.gapA_clip2_norm_y = parameters.get("gapA_clip2_norm_y");
	config.gapA_clip2_norm_z = parameters.get("gapA_clip2_norm_z");	
	
      } // gapA clip control parammeters

    config.do_gapA_slice         = parameters.get_item("do_gapA_slice",{"no","yes"});
    if ( config.do_gapA_slice ==1 )
      {
        config.gapA_slice1_point_x = parameters.get("gapA_slice1_point_x");
	config.gapA_slice1_point_y = parameters.get("gapA_slice1_point_y");
	config.gapA_slice1_point_z = parameters.get("gapA_slice1_point_z");
        config.gapA_slice1_norm_x = parameters.get("gapA_slice1_norm_x");
	config.gapA_slice1_norm_y = parameters.get("gapA_slice1_norm_y");
	config.gapA_slice1_norm_z = parameters.get("gapA_slice1_norm_z");

        // config.gapA_slice2_point_x = parameters.get("gapA_slice2_point_x");
	// config.gapA_slice2_point_y = parameters.get("gapA_slice2_point_y");
	// config.gapA_slice2_point_z = parameters.get("gapA_slice2_point_z");
        // config.gapA_slice2_norm_x = parameters.get("gapA_slice2_norm_x");
	// config.gapA_slice2_norm_y = parameters.get("gapA_slice2_norm_y");
	// config.gapA_slice2_norm_z = parameters.get("gapA_slice2_norm_z");	
	
      } // gapA slice control parammeters

    
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

    config.do_Temperature_clip = parameters.get_item("do_Temperature_clip",{"no", "yes"});
    if (config.do_Temperature_clip ==1)
      {
        config.Temperature_clip_point_x = parameters.get("Temperature_clip_point_x");
	config.Temperature_clip_point_y = parameters.get("Temperature_clip_point_y");
	config.Temperature_clip_point_z = parameters.get("Temperature_clip_point_z");
        config.Temperature_clip_norm_x = parameters.get("Temperature_clip_norm_x");
	config.Temperature_clip_norm_y = parameters.get("Temperature_clip_norm_y");
	config.Temperature_clip_norm_z = parameters.get("Temperature_clip_norm_z");
	config.Temperature_clamp       = parameters.get("Temperature_clamp");
      }

    config.do_Temperature_slice = parameters.get_item("do_Temperature_slice",{"no", "yes"});
    if (config.do_Temperature_slice ==1)
      {
        config.Temperature_slice_point_x = parameters.get("Temperature_slice_point_x");
	config.Temperature_slice_point_y = parameters.get("Temperature_slice_point_y");
	config.Temperature_slice_point_z = parameters.get("Temperature_slice_point_z");
        config.Temperature_slice_norm_x = parameters.get("Temperature_slice_norm_x");
	config.Temperature_slice_norm_y = parameters.get("Temperature_slice_norm_y");
	config.Temperature_slice_norm_z = parameters.get("Temperature_slice_norm_z");
      }

    
    config.do_Temperature_isosurface = parameters.get_item("do_Temperature_isosurface",{"no","yes"});
    if ( config.do_Temperature_isosurface ==1 )
      {
	std::vector<real_t> tmp3 = parameters.get("Temperature_iso_values_vector");
	real_t TcpmK = MP.Tcp_mK(config.Inip);
        for (auto i : tmp3) { config.Temperature_iso_values_vector.push_back(i * TcpmK); }
      } // gapA clip control parammeters

    config.do_phaseMarker_slice = parameters.get_item("do_phaseMarker_slice",{"no","yes"});
    if (config.do_phaseMarker_slice == 1)
      {
        config.pMarker_slice_point_x = parameters.get("pMarker_slice_point_x");
	config.pMarker_slice_point_y = parameters.get("pMarker_slice_point_y");
	config.pMarker_slice_point_z = parameters.get("pMarker_slice_point_z");
        config.pMarker_slice_norm_x = parameters.get("pMarker_slice_norm_x");
	config.pMarker_slice_norm_y = parameters.get("pMarker_slice_norm_y");
	config.pMarker_slice_norm_z = parameters.get("pMarker_slice_norm_z");
      }
    
    // config.do_gapA_3slice       = parameters.get_item("do_gapA_3slice",{"no","yes"});
    // config.do_fe_slice          = parameters.get_item("do_fe_slice",{"no","yes"});
    // config.do_gapA_slice        = parameters.get_item("do_gapA_slice",{"no","yes"});            
    
    config.clamp_bias_gapMin = parameters.get("clamp_bias_gapMin");
    config.clamp_bias_gapMax = parameters.get("clamp_bias_gapMax");
    config.clamp_bias_fed_Min = parameters.get("clamp_bias_fed_Min");
    config.clamp_bias_fed_Max = parameters.get("clamp_bias_fed_Max");
    
    config.camera1_azi = parameters.get("camera1_azi");
    config.camera1_ele = parameters.get("camera1_ele");
    config.camera2_azi = parameters.get("camera2_azi");
    config.camera2_ele = parameters.get("camera2_ele");
    
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


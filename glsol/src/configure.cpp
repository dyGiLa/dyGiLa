#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
//#include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/fft.h"

#include "glsol.hpp"


const std::vector<std::string> glsol::configure(const std::string &fname, int argc, char **argv) {
  
    hila::input parameters(fname);
    hila::out0 << "------------------------------------------------------------" << "\n"
               << "-- dyGiLa 3D p-Wave TDGL HPC Simulation Suites Parameters --" << "\n"
               << "------------------------------------------------------------" << std::endl;
    
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
    config.use_antiQuench         = parameters.get_item("use_antiQuench",{"no", "yes"});
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
    /*
     * gamma2 appling time, started from first moment when blob confCatch() block is entred.
     * This for blob confCatch function noise removment. After this time, T-dependent gamma 
     * is recovery without thermal noise. (in unit of tGL) 
     */ 
    config.extinguish_off_t_count = parameters.get("extinguish_off_t_count");
    
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
	                                                              ,"hotblob"              //8
                                                                      ,"A-n-B"                //9
								      ,"Wiman2016StripeB"});  //10
                                                                       

    hila::out0 << " config.initialCondition is "
	       << config.initialCondition
	       << std::endl;
    
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
	config.Blob_Tc_cutoff = parameters.get_item("Blob_Tc_cutoff",{"no", "yes"});
        config.Ttdb1         = parameters.get("Ttdb1");
        config.Ttdb0         = parameters.get("Ttdb0");
        config.t1            = parameters.get("t1");
       /*********************************************/
       /* >>>>>> hot bloob parameters end here  <<< */    
       /*********************************************/    	
      }

    config.use_CustomerDctxi = parameters.get_item(" use_CustomerDctxi",{"no", "yes"});
    if ( config.use_CustomerDctxi == true ) { hila::out0 << " use Custom Dctxi. " << std::endl; }
    config.Dctxi             = parameters.get("Dctxi");
    
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

    // pStreaming Squezeing Ratio; PSSRatio, the pIO steps skip ratio
    // This control the frequency which pario.pStream() is called;
    // minimum is 1, which means not steps skip, must be unsigned int
    // together with hdf5start/end, one has Precious control on pIO.
    config.PSSRatio = parameters.get("PSSRatio");

    // output_file is the saving path of output file, which offered in congigration file
    // vector for return by names of stream files
    const std::string output_file    = parameters.get("output_file");
    const std::string pVcounter_file = parameters.get("pVcounter_file");

    std::vector<std::string> name_files{output_file, pVcounter_file};

    // xdmf file name, which is provided through config file
    config.xmf2_Amatrix_fname = parameters.get("xmf2_Amatrix_fname");
    config.xmf2_pMarker_fname = parameters.get("xmf2_pMarker_fname");
    config.xmf2_massCurrent_fname = parameters.get("xmf2_massCurrent_fname");
    config.xmf2_spinCurrent_fname = parameters.get("xmf2_spinCurrent_fname");    

    /*----------------------------------------*/
    /* >>>>>>>  boundary conditions  <<<<<<<<<*/
    /*----------------------------------------*/
    config.BCs1 = parameters.get_item("BCs1",{"periodic",       //0
					      "AB",             //1
					      "PairBreaking",   //2
					      "PB_y",           //3
					      "PairB_yz",       //4
					      "BB",             //5
					      "AdGRz",          //6
					      "phaseVortices",  //7
					      "A-n-B"});        //8
    
    config.BCs2 = parameters.get_item("BCs2",{"periodic",       //0
					      "AB",             //1
					      "PairBreaking",   //2
                                              "PB_y",           //3
                                              "PairB_yz",       //4 
                                              "BB",             //5
					      "AdGRz",          //6
                                              "phaseVortices",  //7
                                              "A-n-B"});        //8
    
    config.BCchangec = parameters.get("BCchangec");
    config.use_AdGRz_surfaces = parameters.get_item("use_AdGRz_surfaces",{"no","yes"});
    if ( config.use_AdGRz_surfaces == 1 ) { config.bt = parameters.get("bt_AdGR"); }
    // config.Wn = parameters.get("BoundaryPhaseWindingNO");
           
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
    /* Approx. Gaussian LP filter parameters  */
    /*----------------------------------------*/
    config.useGaussianLP_filter = parameters.get_item("useGaussianLP_filter",{"no","yes"});
    if (config.useGaussianLP_filter == 1)
      {
       config.numIterGLPfilter = parameters.get("numIterGLPfilter");
       // config.GLPfc1 = parameters.get("GLPfc1");
       config.GLPfc2 = parameters.get("GLPfc2");
       /* 3d finite difference coefficients relation */
       config.GLPfc1 = 1. - 6.*config.GLPfc2; 
      }
    
    /*----------------------------------------*/
    /* Parallel IO Engine control parameters  */
    /*----------------------------------------*/
    config.pario_compute_feDensity     = parameters.get_item("pario_compute_feDensity",{"no","yes"});    
    
    config.hdf5_A_matrix_output        = parameters.get_item("hdf5_A_matrix_output",{"no","yes"});    
    // config.hdf5_trA_output             = parameters.get_item("hdf5_trA_output",{"no","yes"});
    // config.hdf5_eigvA_output           = parameters.get_item("hdf5_eigvA_output",{"no","yes"});
    config.hdf5_mass_current_output    = parameters.get_item("hdf5_mass_current_output",{"no","yes"});
    config.hdf5_spin_current_output    = parameters.get_item("hdf5_spin_current_output",{"no","yes"});
    config.hdf5_pMarker_output         = parameters.get_item("hdf5_pMarker_output",{"no","yes"});
    if (
        (config.hdf5_A_matrix_output == 1)
	|| (config.hdf5_mass_current_output == 1)
	|| (config.hdf5_spin_current_output == 1)
	|| (config.hdf5_pMarker_output == 1)
       )
      {
       config.hdf5Ststart = parameters.get("hdf5Ststart");
       config.hdf5Stend   = parameters.get("hdf5Stend");       
      }
    

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

    config.do_phaseMarker_isosurface = parameters.get_item("do_phaseMarker_isosurface",{"no","yes"});
    if (config.do_phaseMarker_isosurface == 1)
      {
	std::vector<real_t> tmp4 = parameters.get("phaseMarker_iso_values_vector");
        for (auto i : tmp4) { config.phaseMarker_iso_values_vector.push_back(i); }
      }

    config.remove_screen_annotations = parameters.get_item("remove_screen_annotations",{"no","yes"});
    config.do_phaseMarker_fieldclip = parameters.get_item("do_phaseMarker_fieldclip",{"no","yes"});
    config.do_phaseMarker_fieldclip_Bphase = parameters.get_item("do_phaseMarker_fieldclip_Bphase",{"no","yes"});
    config.do_phaseMarker_fieldclip_Aphase = parameters.get_item("do_phaseMarker_fieldclip_Aphase",{"no","yes"});    
    
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

    return name_files; //output_file pVfile;
    
} // allocate() function ends here


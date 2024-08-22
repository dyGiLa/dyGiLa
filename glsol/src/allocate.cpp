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


const std::string glsol::allocate(const std::string &fname, int argc, char **argv) {
  
    hila::initialize(argc, argv);
    hila::input parameters(fname);
    
    config.lx = parameters.get("Nx");
    config.ly = parameters.get("Ny");
    config.lz = parameters.get("Nz");
    config.dx = parameters.get("dx");
    config.dtdxRatio = parameters.get("dtdxRatio");
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

    //initialCondition-T
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

    //initialCondition-p
    config.initialConditionp = parameters.get_item("initialConditionp",{"constant"});
    config.Inip	= parameters.get("Inip");
    
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
    if(config.evolveT ==1)
      {
	config.Tevolvetype = parameters.get_item("Tevolvetype",{"heat","wave"});
	config.startdiffT = parameters.get("startdiffT");
	config.diffT = parameters.get("diffT");
      }
    // config.bloob_after = parameters.get_item("bloob_after",{"no","yes"});
    // if(config.bloob_after== 1)
    //   {
    // 	config.theat = parameters.get("theat");
    //   }

    config.useTbath = parameters.get_item("useTbath",{"no","yes"});


    
    /*----------------------------------------*/
    /* Parallel IO Engine control parameters  */
    /*----------------------------------------*/
    config.hdf5_A_matrix_output        = parameters.get_item("hdf5_A_matrix_output",{"no","yes"});
    // config.hdf5_trA_output             = parameters.get_item("hdf5_trA_output",{"no","yes"});
    // config.hdf5_eigvA_output           = parameters.get_item("hdf5_eigvA_output",{"no","yes"});
    config.hdf5_mass_current_output    = parameters.get_item("hdf5_mass_current_output",{"no","yes"});
    config.hdf5_spin_current_output    = parameters.get_item("hdf5_spin_current_output",{"no","yes"});        

    // config.do_gapA_clip         = parameters.get_item("do_gapA_clip",{"no","yes"});
    // config.do_gapA_isosurface   = parameters.get_item("do_gapA_isosurface",{"no","yes"});
    // config.do_gapA_3slice       = parameters.get_item("do_gapA_3slice",{"no","yes"});
    // config.do_fe_slice          = parameters.get_item("do_fe_slice",{"no","yes"});
    // config.do_gapA_slice        = parameters.get_item("do_gapA_slice",{"no","yes"});            
    
    config.clamp_bias_gapMin = parameters.get("clamp_bias_gapMin");
    config.clamp_bias_gapMax = parameters.get("clamp_bias_gapMax");
    config.clamp_fed_Min = parameters.get("clamp_fed_Min");
    config.clamp_fed_Max = parameters.get("clamp_fed_Max");    
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


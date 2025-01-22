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
#include "plumbing/globals.h" 

#include "glsol.hpp"

#include "matep_namespace_utils.hpp" 
//#include "matep.hpp"

#if defined USE_PARIO 
#include "pario.hpp"
#endif

int main(int argc, char **argv) {

    glsol gl;
    
    const std::string output_fname = gl.allocate("sim_params.txt", argc, argv);

    // initialize Temperature field
    gl.initializeT();

    // initialize pressure field
    // gl.initializep();

    // initilize static H-field
    gl.initializeH();

    // initialize OP field
    gl.initialize();

    // initialize blobal matep wrapper
    matep::init_wrapper_mp();    
           
    int stepspos;

    std::initializer_list<int> coordsList {0,0,0};
    const CoordinateVector originpoints(coordsList); 
    
    const unsigned int steps = (gl.config.tEnd - gl.config.tStats) 
                               / (gl.config.dt * gl.config.nOutputs); // number of steps between printing stats
    
    // if (steps == 0)
    //     steps = 1;

    if (gl.config.TDependnetgamma == 0) { gl.config.gamma = gl.config.gamma1; } // initial gamma parameter if gamma is fixed
    
    gl.config.boundaryConditions = gl.config.BCs1;                   // initial bounaryConstions
    
    if (gl.config.positions == 1) {
      stepspos = (gl.config.tEnd - gl.config.tStats) 
                  / (gl.config.dt * gl.config.npositionout);
      if (stepspos == 0)
        stepspos = 1;
    }

    // measurement and stream counter
    unsigned int stat_counter = 0;

    if (hila::myrank() == 0) {
        gl.config.stream.open(output_fname, std::ios::out);
    }

#if defined USE_PARIO
    parIO paraio;
    
    //xml files for MetaData.    
    if (
        (gl.config.hdf5_A_matrix_output        == 1)
	// || (gl.config.hdf5_trA_output          == 1)
	// || (gl.config.hdf5_eigvA_output        == 1)
	|| (gl.config.hdf5_mass_current_output == 1)
	|| (gl.config.hdf5_spin_current_output == 1)
       )
    {
      paraio.xml(gl);
     //hila::synchronize();
    }
    
    paraio.init(gl);

    if (hila::myrank() == 0) {paraio.xdmf(gl);}              
    hila::out0 << "parallel IO enigne starts!" << "\n"
               << std::endl;    
#endif    
    
    /*-------------------------------------------------------------------*/
    /* Dynamic simulation starts after below.                            */
    /*                                                                   */
    /* on gpu the simulation timer is fake, because there's no sync here.*/  
    /* BUt we want to avoid unnecessary sync anyway.                     */
    /*-------------------------------------------------------------------*/      
    static hila::timer run_timer("Simulation time"), meas_timer("Measurements");
    run_timer.start();
    
    while (gl.t < gl.config.tEnd) {
      //gl.config.gamma = (stat_counter < gl.config.gammaoffc) ? gl.config.gamma1 : gl.config.gamma2;
      //hila::out0 << "gl.config.gamma is " << gl.config.gamma << "\n" << std::endl;
      
        if (gl.t >= gl.config.tStats) {
	  
	   if (stat_counter % steps == 0) {

	      meas_timer.start();
	      //gl.write_moduli();
	      gl.write_energies();
	      //gl.write_phases();
	      hila::out0 << "write_energies() call is done "
			 << std::endl;

#if defined USE_PARIO
	      if (
		  (gl.config.hdf5_A_matrix_output == 1)
		  //ToDo: (gl.config.hdf5_mass_current_output == 1) ...
		  //ToDo: (gl.config.hdf5_spin_current_output == 1) ...		  		  
		   && (gl.t >= gl.config.hdf5Ststart && gl.t <= gl.config.hdf5Stend)
		 )
		paraio.pstream(gl);
	      else if (
		       // insitu visualization block, no parallel hd5 stream
                       (gl.config.hdf5_A_matrix_output != 1)
		       && (
                           (gl.config.do_gapA_clip == 1)
			   || (gl.config.do_fed_clip == 1)
			   || (gl.config.do_Temperature_clip == 1)
			   || (gl.config.do_gapA_isosurface == 1)
                          )
                      )
		paraio.pstream(gl);		
	      hila::out0 << "paraio.pstream() call is done " << std::endl;
#endif	            
	      meas_timer.stop();
	   } // streaming block

	   
            /* config.gamma update block */
	    if (
		(gl.config.TDependnetgamma == 1)
		 && (stat_counter < (gl.config.gammaoffc)*steps)
	       )
	      // update gl.config.gamma if T-dependency of gamma is turned on
	      {
		if (gl.config.initialConditionT != 2)
		  {
		   gl.config.gamma = (gl.T.get_element(originpoints) < gl.MP.Tcp_mK(gl.config.Inip))
		     ? gl.MP.gamma_td(gl.config.Inip, gl.T.get_element(originpoints), gl.phaseMarker.get_element(originpoints))
		     : gl.MP.gamma_td(gl.config.Inip, gl.MP.Tcp_mK(gl.config.Inip), gl.phaseMarker.get_element(originpoints));

		  }
		
	      }
	    else if (stat_counter >= (gl.config.gammaoffc)*steps)
      	      // Set gamma to 2nd value after certain momentum no matter what
	      {gl.config.gamma = gl.config.gamma2;}
	    /* config.gamma update block end here */
	    
	    //if (stat_counter == (gl.config.gammaoffc + 3)*steps) {gl.config.gamma = gl.config.gamma1;}
	    if (stat_counter == (gl.config.BCchangec)*steps) {gl.config.boundaryConditions = gl.config.BCs2;}

	    ++stat_counter;

        } //gl.t > gl.config.Stats block

        gl.phaseMarking();        
        hila::out0 << "phaseMarking() call is done. "
    		   << std::endl;
	

	/*******************************************************************/
	/* the following if else blocks are different system synamic updates  
         *  all next_xxx() functions call happen here.
         */
	/*******************************************************************/	
	
	if (
	    (gl.config.initialConditionT == 2)
	    && (gl.config.evolveT == 1)
	   )
	  {
	    //hila::out0 << "just before call next-blob() " << std::endl;
	    gl.next_bath_hotblob_quench_Hfield();
	    hila::out0 << "gl.t is " << gl.t 
		       << ", next_bath_hotblob_quench_Hfield() call " 
		       << ", Tc is " << gl.MP.Tcp_mK(gl.config.Inip)
		       << ", Ttdb0 is " << gl.config.Ttdb0	      
	               << std::endl;
	  } // heterogenous quench, hot blob T-profile 
	else
	  { // homogenous quench block starts from here
	if (
	    ((gl.config.useTbath == 1)
	     && (gl.t >= gl.config.Tbath_start)
	     && (gl.config.evolveT == 0)
	     && (gl.config.gamma.abs() < gl.config.gamma2.abs())) // fixed T-thermal bath
	    ||
	    ((gl.config.useTbath == 1)
	     && (gl.t >= gl.config.Tbath_start)
	     && (gl.config.evolveT == 1)
	     && (gl.t <= gl.config.tThermalizationWaiting)
	     && (gl.config.gamma.abs() < gl.config.gamma2.abs())) // eolved T-thermal bath run, but in fixed T thermalization	    
	   )
	  {
	    gl.next_bath();
	    hila::out0 << " gl.t is " << gl.t << ", gl.config.gamma is " << gl.config.gamma
		       << ", next_bath() call, T in site is " << gl.T.get_element(originpoints)
		       << " Tc is " << gl.MP.Tcp_mK(gl.config.Inip)
		       << std::endl;	    
	  }
	else if (
	         (gl.config.useTbath == 1)
	         && (gl.t >= gl.config.Tbath_start)
		 && (gl.config.evolveT == 1)
		 && (gl.config.Tevolvetype == 2)
		 && (gl.t > gl.config.tThermalizationWaiting)
		 && (gl.config.withHfield != 1)
		 && (gl.config.gamma.abs() < gl.config.gamma2.abs())
                )
	  {
            gl.next_bath_UniT_quench();
	    hila::out0 << " gl.t is " << gl.t << ", gl.config.gamma is " << gl.config.gamma
		       << ", next_bath_UniT_quench() call, T in site is " << gl.T.get_element(originpoints)
		       << ", Tc is " << gl.MP.Tcp_mK(gl.config.Inip)
		       << std::endl;	    

	  }
	else if (
	         (gl.config.useTbath == 1)
	         && (gl.t >= gl.config.Tbath_start)
		 && (gl.config.evolveT == 1)
		 && (gl.config.Tevolvetype == 2)
		 && (gl.t > gl.config.tThermalizationWaiting)
		 && (gl.config.withHfield == 1)
		 && (gl.config.gamma.abs() < gl.config.gamma2.abs())
                )
	  {
            gl.next_bath_UniT_quench_Hfield();
	    hila::out0 << " gl.t is " << gl.t << ", gl.config.gamma is " << gl.config.gamma
		       << ", next_bath_UniT_quench_Hfield() call, T in site is " << gl.T.get_element(originpoints)
	               << ", |H| is " << norm(gl.H.get_element(originpoints))
		       << ", Tc is " << gl.MP.Tcp_mK(gl.config.Inip)
		       << std::endl;	    

	  }
	else
	  {
	    // if gl.config.gamma != gl.config.gamma1, program is doing frozen structure
	    gl.next();
	    hila::out0 << " gl.t is " << gl.t << ", gl.config.gamma is " << gl.config.gamma
		       << ", next() call, T in site is " << gl.T.get_element(originpoints)
		       << ", Tc is " << gl.MP.Tcp_mK(gl.config.Inip)
		       << std::endl;	    
	    
	  }

	  } // homogenous quench block
	    
		
    } // gl.t evolves while loop ends here
    run_timer.stop();

#if defined USE_PARIO
    paraio.shutdown();
    hila::out0 << "parallel IO engine shutdown! " << "\n"
               << std::endl;    
#endif    

    if (hila::myrank() == 0) {
        gl.config.stream.close();
    }

    hila::finishrun();
    return 0;
}

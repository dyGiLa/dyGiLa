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
#include "he3sim.hpp"
#include "pario.hpp"
#include "matep.hpp"

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

int main(int argc, char **argv) {
    he3sim sim;
    const std::string output_fname = sim.allocate("sim_params.txt", argc, argv);
    sim.update_params();
    sim.initialize();
    int stepspos;
    
    int steps = (sim.config.tEnd - sim.config.tStats)
                 / (sim.config.dt * sim.config.nOutputs); // number of steps between out stream.

    if (steps == 0)
        steps = 1;

    sim.config.gamma              = sim.config.gamma1;                 // initial gamma parameter
    sim.config.boundaryConditions = sim.config.BCs1;                   // initial bounaryConstions
    /*if (sim.config.positions == 1)
      {
	stepspos =  (sim.config.tEnd - sim.config.tStats)
	             /(sim.config.dt * sim.config.npositionout);
	if (stepspos == 0)
          stepspos = 1;
	  }*/
	
    int stat_counter = 0;

    if (hila::myrank() == 0) {
        sim.config.stream.open(output_fname, std::ios::out);
    }

    // xmf2 file output for paraview.    
    if (
        (sim.config.hdf5_A_matrix_output        == 1)
	|| (sim.config.hdf5_trA_output          == 1)
	|| (sim.config.hdf5_eigvA_output        == 1)
	|| (sim.config.hdf5_mass_current_output == 1)
	|| (sim.config.hdf5_spin_current_output == 1)
       ){
     sim.insitu_hdf5xdmf();
     //hila::synchronize();    
     }

    if (hila::myrank() == 0) sim.write_xdmf();    
    
#if defined USE_ASCENT
    hila::out0 << "using ascent" << "\n\n\n";
    sim.insitu_initialize();
#endif    
    
    /*-------------------------------------------------------------------*/
    /* Dynamic simulation starts after below.                            */
    /*                                                                   */
    /* on gpu the simulation timer is fake, because there's no sync here.*/  
    /* BUt we want to avoid unnecessary sync anyway.                     */
    /*-------------------------------------------------------------------*/
    static hila::timer run_timer("Simulation time"), meas_timer("Measurements");
    run_timer.start();
    

    while (sim.t < sim.config.tEnd) {
      //sim.config.gamma = (stat_counter < sim.config.gammaoffc) ? sim.config.gamma1 : sim.config.gamma2;
      //hila::out0 << " sim_config.gamma is " << sim.config.gamma << "\n" << std::endl;
      
        if (sim.t >= sim.config.tStats) {
	  
            if (stat_counter % steps == 0) {

	      meas_timer.start();
	      //sim.write_moduli();
	      sim.write_energies();
	      //sim.write_phases();

#if defined USE_ASCENT
              sim.insitu_execute();
#endif	      
	      meas_timer.stop();

            }
	    if (sim.config.positions == 1)
	      {
		if (stat_counter % stepspos == 0)
		  {
		  //sim.write_A_matrix_positions();
                   sim.write_positions();
		  }
	      }
	    if (stat_counter == (sim.config.gammaoffc)*steps) {sim.config.gamma = sim.config.gamma2;}
	    //if (stat_counter == (sim.config.gammaoffc + 3)*steps) {sim.config.gamma = sim.config.gamma1;}
	    if (stat_counter == (sim.config.BCchangec)*steps) {sim.config.boundaryConditions = sim.config.BCs2;}
            stat_counter++;
        } //sim.t > sim.config.Stats block	

	sim.update_params();
	
	if (sim.config.useTbath == 1)
	  {
	    sim.next_bath();
	  }
	else
	  {
	    sim.next();
	  }
    } // wile loop sim.t ends here
    run_timer.stop();

#if defined USE_ASCENT
    sim.insitu_close();
    hila::out0 << "ascent closed." << "\n\n\n";    
#endif    

    if (hila::myrank() == 0) {
        sim.config.stream.close();
    }

    hila::finishrun();
    return 0;
}

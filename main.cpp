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
#include "pario.hpp"

int main(int argc, char **argv) {

    glsol gl;
    parIO paraio;
    
    const std::string output_fname = gl.allocate("sim_params.txt", argc, argv);
    gl.update_params();
    gl.initialize();
    
    int stepspos;
    
    int steps = (gl.config.tEnd - gl.config.tStats)
                 / (gl.config.dt * gl.config.nOutputs); // number of steps between out stream.

    if (steps == 0)
        steps = 1;

    gl.config.gamma              = gl.config.gamma1;                 // initial gamma parameter
    gl.config.boundaryConditions = gl.config.BCs1;                   // initial bounaryConstions
    /*if (sim.config.positions == 1)
      {
	stepspos =  (sim.config.tEnd - sim.config.tStats)
	             /(sim.config.dt * sim.config.npositionout);
	if (stepspos == 0)
          stepspos = 1;
	  }*/
	
    int stat_counter = 0;

    if (hila::myrank() == 0) {
        gl.config.stream.open(output_fname, std::ios::out);
    }

    // xmf2 file output for paraview.    
    // if (
    //     (gl.config.hdf5_A_matrix_output        == 1)
    // 	|| (gl.config.hdf5_trA_output          == 1)
    // 	|| (gl.config.hdf5_eigvA_output        == 1)
    // 	|| (gl.config.hdf5_mass_current_output == 1)
    // 	|| (gl.config.hdf5_spin_current_output == 1)
    //    )
     {
      paraio.xdmf(gl);
     //hila::synchronize();    
     }

    if (hila::myrank() == 0) gl.write_xdmf();    
    
#if defined USE_PARIO
    hila::out0 << "parallel IO enigne starts!" << "\n\n\n";
    paraio.init(gl);
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
      //sim.config.gamma = (stat_counter < sim.config.gammaoffc) ? sim.config.gamma1 : sim.config.gamma2;
      //hila::out0 << " sim_config.gamma is " << sim.config.gamma << "\n" << std::endl;
      
        if (gl.t >= gl.config.tStats) {
	  
            if (stat_counter % steps == 0) {

	      meas_timer.start();
	      //sim.write_moduli();
	      gl.write_energies();
	      //sim.write_phases();

#if defined USE_PARIO
              paraio.pstream(gl);
#endif	      
	      meas_timer.stop();

            }
	    // if (gl.config.positions == 1)
	    //   {
	    // 	if (stat_counter % stepspos == 0)
	    // 	  {
	    // 	  //sim.write_A_matrix_positions();
            //        gl.write_positions();
	    // 	  }
	    //   }
	    if (stat_counter == (gl.config.gammaoffc)*steps) {gl.config.gamma = gl.config.gamma2;}
	    //if (stat_counter == (sim.config.gammaoffc + 3)*steps) {sim.config.gamma = sim.config.gamma1;}
	    if (stat_counter == (gl.config.BCchangec)*steps) {gl.config.boundaryConditions = gl.config.BCs2;}
            stat_counter++;

	} //sim.t > sim.config.Stats block	

	gl.update_params();
	
	if (gl.config.useTbath == 1)
	  {
	    gl.next_bath();
	  }
	else
	  {
	    gl.next();
	  }
    } // wile loop sim.t ends here
    run_timer.stop();

#if defined USE_PARIO
    paraio.shutdown();
    hila::out0 << "parallel IO engine shutdown! " << "\n\n\n";    
#endif    

    if (hila::myrank() == 0) {
        gl.config.stream.close();
    }

    hila::finishrun();
    return 0;
}

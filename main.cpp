#define USE_PARIO 
#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/fft.h"
//#include "plumbing/globals.h" 

#include "glsol.hpp"
//#include "matep_namespace_utils.hpp"
#include "orch.hpp"

#if defined USE_PARIO 
#include "pario.hpp"
#endif

int main(int argc, char **argv) {

    glsol gl;      
       
    const std::vector<std::string> name_files = gl.allocate("sim_params.txt", argc, argv);
    // initialize blobal matep wrapper
    matep::init_wrapper_mp();           
    // initialize Temperature field
    gl.initializeT();
    // initialize pressure field
    // gl.initializep();
    // initilize static H-field
    gl.initializeH();
    // initialize OP field
    gl.initialize();   
           
    std::initializer_list<int> coordsList {0,0,0};
    const CoordinateVector originpoints(coordsList); 

    // number of steps between reduction streaming
    const unsigned int steps = (gl.config.tEnd - gl.config.tStats) 
                               /(gl.config.dt * gl.config.nOutputs); 
    
    // initial gamma parameter if gamma is fixed
    if (gl.config.TDependnetgamma == 0) { gl.config.gamma = gl.config.gamma1; } 

    // initial bounaryConstions
    gl.config.boundaryConditions = gl.config.BCs1;
    if (hila::myrank() == 0) { gl.fstreams_open(name_files); }    
        
#if defined USE_PARIO
    parIO paraio;    
    //xml files for MetaData.    
    orch::writeHDF5_xmls(gl, paraio);    
    paraio.init(gl);

    if (hila::myrank() == 0) paraio.xdmf(gl);              
    hila::out0 << "parallel IO enigne starts!" << std::endl;    
#endif    
    
    /*-------------------------------------------------------------------*/
    /* Dynamic simulation starts after below.                            */
    /*                                                                   */
    /* on gpu the simulation timer is fake, because there's no sync here.*/  
    /* But we want to avoid unnecessary sync anyway.                     */
    /*-------------------------------------------------------------------*/      
    static hila::timer run_timer("Simulation time"), meas_timer("Measurements");
    run_timer.start();
    // measurement and stream counter
    unsigned int stat_counter = 0;
        
    while (gl.t < gl.config.tEnd) {
      
        if (gl.t > gl.config.tStats) {
	  
	   if (stat_counter % steps == 0) {
	      meas_timer.start();
	      orch::pStreaming(gl, paraio, stat_counter);
	      meas_timer.stop();
	   } // streaming block

	   // gamma handling
	   orch::gammaEvolve(gl, stat_counter, steps, originpoints);
	   
	   if (stat_counter == (gl.config.BCchangec)*steps)
	     { gl.config.boundaryConditions = gl.config.BCs2; }

	   ++stat_counter;

        } //gl.t > gl.config.Stats block

	if (gl.config.TDependnetgamma == true) {
	    // do every-dt phase-Marking when gamma is T-dependent
            gl.phaseMarking();        
            hila::out0 << "gl.t is " << gl.t
		       << ", phaseMarking() call is done. "
		       << std::endl;
	  }

        // t-evolve call
        orch::nextBlocks(gl, stat_counter, steps, originpoints);
         	    		
    } // gl.t evolves while loop ends here
    run_timer.stop();

#if defined USE_PARIO
    paraio.shutdown();
    hila::out0 << "parallel IO engine shutdown! " << std::endl;    
#endif    

    if (hila::myrank() == 0) { gl.fstreams_close(); }

    hila::finishrun();
    return 0;
}

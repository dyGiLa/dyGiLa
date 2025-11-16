#define USE_PMD_GAMMA
#define USE_PARIO 
#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>

#include "plumbing/hila.h"
// #include "plumbing/fft.h"
#include "plumbing/globals.h" 

#include "glsol.hpp"
//#include "matep_namespace_utils.hpp"
#include "dyGiLa.hpp"
#include "dyGiLa_utils.hpp"

#if defined USE_PARIO 
#include "pario.hpp"
#endif

int main(int argc, char **argv) {

    /*----------------------------*/
    /*--- dyGiLa initialization --*/
    /*----------------------------*/
    glsol gl;    
    auto dyGiLa = dyGiLa::dyGiLaInit(gl, argc, argv);
    const std::vector<std::string> name_files = std::get<0>(dyGiLa);
    const CoordinateVector originpoints = *std::get<1>(dyGiLa);
    const unsigned int steps = std::get<2>(dyGiLa);        
  
    // initial gamma parameter if gamma is constant
    if (gl.config.TDependnetgamma == 0) { gl.config.gamma = gl.config.gamma1; } 

    // initial bounaryConstions
    gl.config.boundaryConditions = gl.config.BCs1;
    if (hila::myrank() == 0) { gl.fstreams_open(name_files); }    
        
#if defined USE_PARIO
    parIO paraio;    
    //xml files for MetaData.    
    dyGiLa::writeHDF5_xmls(gl, paraio);    
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
#ifdef USE_PMD_GAMMA
	  if (gl.config.TDependnetgamma == true)
	    { dyGiLa::utils::phaseMarking(gl, stat_counter, steps); }
#endif	
	   if (stat_counter % steps == 0) {
	      meas_timer.start();
#ifndef USE_PMD_GAMMA
	      if (gl.config.TDependnetgamma == true)
		{ dyGiLa::utils::phaseMarking(gl, stat_counter, steps); }
#endif	
	      dyGiLa::pStreaming(gl, paraio, stat_counter, steps);
	      meas_timer.stop();
	   } // streaming block

	   // gamma handling
	   dyGiLa::gammaEvolve(gl, stat_counter, steps, originpoints);	   
	   if (stat_counter == (gl.config.BCchangec)*steps)
	     { gl.config.boundaryConditions = gl.config.BCs2; }

	   ++stat_counter;
        } //gl.t > gl.config.Stats block

        // t-evolve call
        dyGiLa::nextBlocks(gl, stat_counter, steps, originpoints);
         	    		
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

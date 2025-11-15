//#define USE_PARIO 
#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/globals.h" 

#include "glsol.hpp"
//#include "matep_namespace_utils.hpp"
#include "dyGiLa.hpp"
#include "dyGiLa_utils.hpp"

// #if defined USE_PARIO 
// #include "pario.hpp"
// #endif

namespace dyGiLa {
  
void
utils::homogenousQuench(glsol &gl, unsigned int &stat_counter, const unsigned int &steps, const CoordinateVector &originpoints, const unsigned int &modSteps, const unsigned int &modPSSR)
{  
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
     && (gl.config.gamma.abs() < gl.config.gamma2.abs())) // eolved T-thermal bath run, but in fixed T thermalization w/o AdGR boundary
   )
  {
    if (gl.config.use_AdGRz_surfaces != 1)
      {
       gl.next_bath();
       if (
           (modSteps == 0)
           && (modPSSR == 0)	   
          )
	 {
          hila::out0 << " gl.t is " << gl.t << ", gl.config.gamma is " << gl.config.gamma
	             << ", next_bath() call, T in site is " << gl.T.get_element(originpoints)
	             << " Tc is " << gl.MP.Tcp_mK(gl.config.Inip)
	             << std::endl;	    
	 }       	          	          	        
      }
    else if (gl.config.use_AdGRz_surfaces == 1)
      {
       gl.next_AdGRz_bath();
       if (
           (modSteps == 0)
           && (modPSSR == 0)	   
	  )
	 {
           hila::out0 << " gl.t is " << gl.t << ", gl.config.gamma is " << gl.config.gamma
	              << ", next_AdGRz_bath() call, T in site is " << gl.T.get_element(originpoints)
	              << " Tc is " << gl.MP.Tcp_mK(gl.config.Inip)
	              << std::endl;	      
	 }            
      }
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
    if (
        (modSteps == 0)
        && (modPSSR == 0)	   	
       )
      {
       hila::out0 << " gl.t is " << gl.t << ", gl.config.gamma is " << gl.config.gamma
	          << ", next_bath_UniT_quench() call, T in site is " << gl.T.get_element(originpoints)
	          << ", Tc is " << gl.MP.Tcp_mK(gl.config.Inip)
	          << std::endl;	    
      }    
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
    if (gl.config.use_AdGRz_surfaces != 1)
      {
       gl.next_bath_UniT_quench_Hfield();
       if (
           (modSteps == 0)
           && (modPSSR == 0)	   	   
          )
	 {
       	  hila::out0 << " gl.t is " << gl.t << ", gl.config.gamma is "
		     << gl.config.gamma
      	             << ", next_bath_UniT_quench_Hfield() call, T in site is "
		     << gl.T.get_element(originpoints)
                     << ", |H| is " << norm(gl.H.get_element(originpoints))
	             << ", Tc is " << gl.MP.Tcp_mK(gl.config.Inip)
	             << std::endl;	    
	 }       
      }
    else if (gl.config.use_AdGRz_surfaces == 1)
      {
        gl.next_bath_UniT_quench_AdGRz_Hfield();
        if (
            (modSteps == 0)
            && (modPSSR == 0)	   	    
           )
	  {
       	   hila::out0 << " gl.t is " << gl.t << ", gl.config.gamma is "
		      << gl.config.gamma
      	              << ", next_bath_UniT_quench_AdGR_Hfield() call, T in site is "
		      << gl.T.get_element(originpoints)
                      << ", |H| is " << norm(gl.H.get_element(originpoints))
	              << ", Tc is " << gl.MP.Tcp_mK(gl.config.Inip)
	              << std::endl;	    
	  }
	
      }

  }
else if (
         (gl.config.useTbath == 1)
         && (gl.t >= gl.config.Tbath_start)
	 && (gl.config.evolveT == 1)
	 && (gl.config.Tevolvetype == 2)
	 && (gl.t > gl.config.tThermalizationWaiting)
	 && (gl.config.withHfield == 1)
	 && (gl.config.gamma.abs() >= gl.config.gamma2.abs())
        )
  { // configuration catch block 
    if (gl.config.use_AdGRz_surfaces == false)
      { /* free UniT quench Hfield */ }
    else if (gl.config.use_AdGRz_surfaces == true)
      {
       ++gl.extinguish_t; //estinguish time count, in step of dt
       gl.next_bath_UniT_quench_AdGRz_Hfield_confCatch();
       if (
           (modSteps == 0)
	   && (modPSSR == 0)	   
          )
         {// squeze IO a little bit
       	   hila::out0 << " gl.t is " << gl.t << ", gl.config.gamma is "
		      << gl.config.gamma
      	              << ", next_bath_UniT_quench_AdGRz_Hfield_confCatch() call, T in site is "
		      << gl.T.get_element(originpoints)
                      << ", |H| is " << norm(gl.H.get_element(originpoints))
	              << ", Tc is " << gl.MP.Tcp_mK(gl.config.Inip)
	              << std::endl;	    	   
         }       
      }
  } // UniT quench configuration catch block ends here 
else
  {
    // if gl.config.gamma != gl.config.gamma1, program is doing frozen structure
    gl.next();
    if (
        (modSteps == 0)
        && (modPSSR == 0)	   	
       )
      {
       hila::out0 << " gl.t is " << gl.t
		  << ", gl.config.gamma is "
		  << gl.config.gamma
	          << ", next() call, T in site is "
		  << gl.T.get_element(originpoints)
	          << ", Tc is " << gl.MP.Tcp_mK(gl.config.Inip)
	          << std::endl;	    
      }   	    
  }

} // homogenous quench func ends here

} // orch namespace ends here

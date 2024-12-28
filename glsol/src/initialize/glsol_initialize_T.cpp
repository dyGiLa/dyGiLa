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


void glsol::initializeT() {

  switch (config.initialConditionT) {

  case 0: {

    // constant uniform temperature
    onsites(ALL) {
      T[X] = config.IniT;
      dT[X] = 0.0;
    }

    hila::out0 << "Uniform T field is initialized. "
               << std::endl;;
    
    break;
  }
   
  case 1: {

    // case of sine configuration
    onsites(ALL){
      
      auto xcoord = X.coordinate(e_x);

      T[X] = config.IniT + config.ampT*sin(2.0*M_PI*xcoord/config.lx);
      dT[X] = 0.0;
    }

     hila::out0 << "Sine T \n";
    
    break;
  }

  case 2: {

    // case of hot bloob
    /* initialize spherical symmetric hot bloob
     * at the momentum when radius of Tc fronter achieves maximum  
     */
    onsites(ALL){
      /* hila's coordinate index is counted from zero at corner,
       * so for a blob at the center of box, you need coordinate transformation
       */
      auto x = X.coordinate(e_x) - config.lx/2.0;
      auto y = X.coordinate(e_y) - config.ly/2.0;
      auto z = X.coordinate(e_z) - config.lz/2.0;

      real_t tm  = MP.t_TcMax_blob(config.Inip, config.Ttdb1, config.Ttdb0, config.t1);
      real_t Tcp_mK = MP.Tcp_mK(config.Inip);
      
      auto r2 = x*x + y*y + z*z;
      
      T[X] = ((config.Ttdb1 - config.Ttdb0) * Tcp_mK
	      * std::pow(config.t1/tm, 3./2.)
	      * exp(-r2/(4. * MP.Dd(config.Inip) * tm)))
	     + config.Ttdb0 * Tcp_mK;

    }

    hila::out0 << "tm is " << MP.t_TcMax_blob(config.Inip, config.Ttdb1, config.Ttdb0, config.t1) * MP.tGL(config.Inip) << "s"
               << "\n"
               << "tGL is " << MP.tGL(config.Inip) << "s"
               << "\n"
               << "t1 is "  << config.t1 * MP.tGL(config.Inip) << "s"
               << "\n"
               << "Dd is "  << MP.Dd(config.Inip) * (MP.xi0GLp(config.Inip) * MP.xi0GLp(config.Inip))/MP.tGL(config.Inip) * (1e6) << "mu-m^2.mu-s^-1"
	       << "\n"
               << "xi0GLp is " << MP.xi0GLp(config.Inip)
               << "\n"
               << "Hot Bloob with radius "
               << MP.r_TcMax_blob(config.Inip, config.Ttdb1, config.Ttdb0, config.t1) * (MP.xi0GLp(config.Inip)) * (1e6) << "mu-m"
               << " in initialized. Tc frontier vasnished at "
               << MP.t_TcVanish_blob(config.Inip, config.Ttdb1, config.Ttdb0, config.t1)      
               << std::endl;
    
    break;
  }  // hot bloob block ends here

 } // switch config.initialConfiguraionT block end
}


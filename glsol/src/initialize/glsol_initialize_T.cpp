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

#ifdef T_FIELD
  switch (config.initialConditionT) {

  case 0: {

    // constant uniform temperature
    onsites(ALL) {
      T[X] = config.IniT;
      dT[X] = 0.0;
    }

    hila::out0 << "Constant T to T="<<config.IniT<<"mK \n";
    
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

    // case of hotbloob 
    onsites(ALL){
      auto xcoord = X.coordinate(e_x);
      auto ycoord = X.coordinate(e_y);
      auto zcoord = X.coordinate(e_z);

      real_t expx=exp(-0.5*pow(xcoord-config.lx/2.0,2.0)/pow(config.sigTx,2.0));
      real_t expy=exp(-0.5*pow(ycoord-config.ly/2.0,2.0)/pow(config.sigTy,2.0));
      real_t expz=exp(-0.5*pow(zcoord-config.lz/2.0,2.0)/pow(config.sigTz,2.0));
      
      T[X] = config.IniT + config.ampT*expx*expy*expz;
      dT[X] = 0.0;
    }

     hila::out0 << "Hot Spot \n";
    
    break;
  } 
  } // switch config.initialConfiguraionT block end
#else
  T=config.IniT;
#endif  
}


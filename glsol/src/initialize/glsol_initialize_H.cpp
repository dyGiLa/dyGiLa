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


void glsol::initializeH() {

  switch (config.initialConditionH) {

  case 0: {

    onsites(ALL) {
      foralldir(al) {
	H[X].e(al) = config.InitH.e(al);
      }

    }

    hila::out0 << "Constant H-field intialized" << std::endl;

     break;
  }
  } // switch config.initialConfiguraionp block end
}



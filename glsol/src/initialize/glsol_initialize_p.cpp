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


void glsol::initializep() {

  switch (config.initialConditionp) {

  case 0: {

    onsites(ALL) {

      p[X] = config.Inip;

    }

     hila::out0 << "Constant p \n";

     break;
  }
  } // switch config.initialConfiguraionp block end
}



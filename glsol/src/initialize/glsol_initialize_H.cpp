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

  H[0] = config.IniHx;
  H[1] = config.IniHy;
  H[2] = config.IniHz;

  hila::out0 << "Magnetic field initialized to: H="<<H[0]<<","<<H[1]<<","<<H[2]<<"\n";
}



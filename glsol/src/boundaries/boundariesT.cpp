#define USE_MPI
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>

#include "plumbing/hila.h"
#include "plumbing/fft.h"

#include "glsol.hpp"
#include "matep.hpp"

real_t glsol::periodic_T(real_t T1) {

  return T1;
}

real_t glsol::fixedT() {

  real_t Tbc;

  Tbc=config.T_boundary;

  return Tbc;

}

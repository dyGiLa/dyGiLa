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
#include "matep.hpp"

void glsol::point_params(real_t T, real_t p, real_t beta[6]) {
  
  beta[0] = MP.alpha_td(p, T);
  beta[1] = MP.beta1_td(p, T);
  beta[2] = MP.beta2_td(p, T);
  beta[3] = MP.beta3_td(p, T);
  beta[4] = MP.beta4_td(p, T);
  beta[5] = MP.beta5_td(p, T);
}
  

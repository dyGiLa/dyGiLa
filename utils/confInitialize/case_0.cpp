#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>

#include "plumbing/hila.h"
// #include "plumbing/fft.h"

#include "glsol.hpp"
#include "matep.hpp"

void glsol::case_0() {
   
    pi = 0;
    deltaPi =0;
    djAaj = 0;
    real_t gap = MP.gap_B_td(config.Inip, config.IniT);
    onsites(ALL) {                     
      A[X] = hila::gaussrand();
      A[X] = gap * A[X]/A[X].norm();   
    }

    hila::out0 << "Components randomly created! " << std::endl;

} // case_0() call end here


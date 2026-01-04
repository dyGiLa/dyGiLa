#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/fft.h"

#include "glsol.hpp"
#include "matep.hpp"


void glsol::case_2() {
       
    pi = 0;
    deltaPi =0;
    djAaj = 0;    
    phaseMarker = 0.F;
    
    onsites(ALL) { A[X] = sqrt(0.1) * hila::gaussrand(); }

    hila::out0 << " normal-phase-real-1 created \n";

} // case_2() call end here


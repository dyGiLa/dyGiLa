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


void glsol::case_3() {
   
    pi = 0;
    deltaPi = 0;
    djAaj = 0;    
    phaseMarker = 0.F;
    
    onsites(ALL) {
      /*foralldir(al) foralldir(i){
        A[X].e(al,i).real().gaussian_random();
	}*/
      A[X].gaussian_random();
      A[X] = A[X].real();
    }
    hila::out0 << " normal-phase-real-2 created! " << std::endl;

} // case_3() call end here


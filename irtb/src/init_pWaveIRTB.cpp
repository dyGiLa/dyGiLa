#include <iostream>
#include <cstddef>
#include <cmath>
#include <vector>

#include "irtb.hpp"
#include "plumbing/globals.h"

//hila::global<irtb::irtb_matrix> pWaveIRTB;    

namespace irtb {

  void init_pWaveIRTB() {

    irtb_matrix IRTB_Matrix;

    // real_t c1_ARR[18] = {-0.0098, -0.0127, -0.0155, -0.0181, -0.0207, -0.0231, -0.0254, -0.0275, -0.0295, -0.0314, -0.0330, -0.0345, -0.0358, -0.0370, -0.0381, -0.0391, -0.0402, -0.0413};
    // for (unsigned int i = 0; i<18; ++i) { mp_consts.c1_arr[i] = c1_ARR[i]; }

    // ...
    
    // assign to global wrapper
    pWaveIRTB = IRTB_Matrix;
    
    
  } // init_pWaveIRTB() functon block ends here
  
} // namespace irtb block ends here

#ifndef IRTB_HPP
#define IRTB_HPP

#include <cstddef>
#include <cmath>
#include <vector>

#include "plumbing/hila.h"
#include "plumbing/globals.h"

namespace irtb {

  //extern template class hila::global<irtb_consts>;
  hila::global<irtb::irtb_matrix> pWaveIRTB;

  struct irtb_matrix {

    Matrix<3,3,Complex<float>> T00;
    Matrix<3,3,Complex<float>> T10;
    Matrix<3,3,Complex<float>> T1p1;
    Matrix<3,3,Complex<float>> T1n1;
    Matrix<3,3,Complex<float>> T20;    
    Matrix<3,3,Complex<float>> T2p1;
    Matrix<3,3,Complex<float>> T2n1;
    Matrix<3,3,Complex<float>> T2p2;
    Matrix<3,3,Complex<float>> T2n2;        

  };

  // function initializing hila::global<irtb::irtb_matrix> pWaveIRTB
  void init_pWaveIRTB();
    
}  // irtb namespace block ends here
#endif

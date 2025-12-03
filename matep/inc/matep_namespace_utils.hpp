/*
 * 
 * This is the *.hpp header file of the Strong Coupling Correction Object (SCCO) 
 * or Material Parameters Object (Matep).
 * 
 * Member functions are declared at here, and then defined at Matep.cpp.
 *
 * author: Quang. Zhang (timohyva@github)
 *
 */ 

#ifndef MATEP_NAMESPACE_UTILS_HPP
#define MATEP_NAMESPACE_UTILS_HPP

#include <iostream>
#include <cstddef>
#include <cmath>
#include <vector>

//#include "matep_consts.hpp"

#include "plumbing/hila.h"
#include "plumbing/globals.h"

namespace matep {

  using real_t = float; // or double

  //extern template class hila::global<matep_consts>;

  struct matep_consts {

  double Kelvin,
         J,
         s,
         m,
         kg,
         pi,
         E,
         p_pcp;

  // physical constants for he3
  double u,
         m3,
         nm,
         hbar,
         kb,
         zeta3,
         c_betai,
         gammahbar,
         tau0N;

  double c1_arr[18];
  double c2_arr[18];
  double c3_arr[18];
  double c4_arr[18];
  double c5_arr[18];

  // ***************************************************
  double Fa0_arr[18];// Landau parameter  
  double Tc_arr[18]; // mK;
  double Ms_arr[18]; // in unit of helium-3 atom;
  double VF_arr[18]; // fermi velosity, m.s^-1;
  double XI0_arr[18];

  };

  // function initializing hila::global<matep::matep_consts> wrapper_mp
  // wrapper_mp is visible for this function, following the inner-outer scope rules
  // wrapper_mp is decleared in outer scope of matep.cpp
  void init_wrapper_mp();
    
}  // matep namespace block ends here
#endif

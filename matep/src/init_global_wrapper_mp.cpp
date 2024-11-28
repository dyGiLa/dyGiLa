/*
 * This is the *.cpp file of the Strong Coupling Correction Object (SCCO)
 * or Material Parameters Object (Matep).
 * 
 * Member functions are declared at Matep.hpp, and then defined at here.
 *
 * author: Quang. Zhang (timohyva@github)
 *
 */ 


#include <iostream>
#include <cstddef>
#include <cmath>
#include <vector>

#include "matep_namespace_utils.hpp"
//#include "matep.hpp"
//#include "wrapper_mp_instantiation_def.hpp" 
//#include "plumbing/hila.h"


//#include "matep_consts.hpp"
#include "plumbing/globals.h"

// extern declearation of wrapper_mp, which is defined in matep.cpp file.
// This proctice prevent multiple defiation of same object wrapper, while still keeps access of single instance of wrapper_mp
extern hila::global<matep::matep_consts> wrapper_mp;    

namespace matep {

  void init_wrapper_mp() {

    matep_consts mp_consts;

    mp_consts.Kelvin = 1.f;
    mp_consts.J = 1.f;
    mp_consts.s = 1.f;
    mp_consts.m = 1.f;
    mp_consts.kg = 1.f;
    mp_consts.pi = 3.14159265358979323846264338328f;
    mp_consts.p_pcp = 21.22f;

    mp_consts.u = 1.66053906660f*(10.0e-27)*mp_consts.kg;
    mp_consts.m3 = 3.016293f*mp_consts.u;
    mp_consts.nm = (10.0e-9)*mp_consts.m;
    mp_consts.hbar = 1.054571817f*(10.0e-34)*mp_consts.J*mp_consts.s;
    mp_consts.kb = 1.380649*(10.0e-23)*mp_consts.J*1.0f;
    mp_consts.zeta3 = 1.2020569031595942;
    mp_consts.c_betai = (7.0f*mp_consts.zeta3)/(80.0f*mp_consts.pi*mp_consts.pi);

    real_t c1_ARR[18] = {-0.0098, -0.0127, -0.0155, -0.0181, -0.0207, -0.0231, -0.0254, -0.0275, -0.0295, -0.0314, -0.0330, -0.0345, -0.0358, -0.0370, -0.0381, -0.0391, -0.0402, -0.0413};
    for (unsigned int i = 0; i<18; ++i) { mp_consts.c1_arr[i] = c1_ARR[i]; }

    real_t c2_ARR[18] = {-0.0419, -0.0490, -0.0562, -0.0636, -0.0711, -0.0786, -0.0861, -0.0936, -0.1011, -0.1086, -0.1160, -0.1233, -0.1306, -0.1378, -0.1448, -0.1517, -0.1583, -0.1645};
    for (unsigned int i = 0; i<18; ++i) { mp_consts.c2_arr[i] = c2_ARR[i]; }

    real_t c3_ARR[18] = {-0.0132, -0.0161, -0.0184, -0.0202, -0.0216, -0.0226, -0.0233, -0.0239, -0.0243, -0.0247, -0.0249, -0.0252, -0.0255, -0.0258, -0.0262, -0.0265, -0.0267, -0.0268};
    for (unsigned int i = 0; i<18; ++i) { mp_consts.c3_arr[i] = c3_ARR[i]; }

    real_t c4_ARR[18] = {-0.0047, -0.0276, -0.0514, -0.0760, -0.1010, -0.1260, -0.1508, -0.1751, -0.1985, -0.2208, -0.2419, -0.2614, -0.2795, -0.2961, -0.3114, -0.3255, -0.3388, -0.3518};
    for (unsigned int i = 0; i<18; ++i) { mp_consts.c4_arr[i] = c4_ARR[i]; }

    real_t c5_ARR[18] = {-0.0899, -0.1277, -0.1602, -0.1880, -0.2119, -0.2324, -0.2503, -0.2660, -0.2801, -0.2930, -0.3051, -0.3167, -0.3280, -0.3392, -0.3502, -0.3611, -0.3717, -0.3815};
    for (unsigned int i = 0; i<18; ++i) { mp_consts.c5_arr[i] = c5_ARR[i]; }


    real_t Tc_ARR[18] = {0.929, 1.181, 1.388, 1.560, 1.705, 1.828, 1.934, 2.026, 2.106, 2.177, 2.239, 2.293, 2.339, 2.378, 2.411, 2.438, 2.463, 2.486}; // mK
    for (unsigned int i = 0; i<18; ++i) { mp_consts.Tc_arr[i] = Tc_ARR[i]; }

    real_t Ms_ARR[18] = {2.80, 3.05, 3.27, 3.48, 3.68, 3.86, 4.03, 4.20, 4.37, 4.53, 4.70, 4.86, 5.02, 5.18, 5.34, 5.50, 5.66, 5.82}; // in unit of helium-3 atom;
    for (unsigned int i = 0; i<18; ++i) { mp_consts.Ms_arr[i] = Ms_ARR[i]; }

    real_t VF_ARR[18] = {59.03, 55.41, 52.36, 49.77, 47.56, 45.66, 44.00, 42.51, 41.17, 39.92, 38.74, 37.61, 36.53, 35.50, 34.53, 33.63, 32.85, 32.23}; // fermi velosity, m.s^-1;
    for (unsigned int i = 0; i<18; ++i) { mp_consts.VF_arr[i] = VF_ARR[i]; }

    real_t XI0_ARR[18] = {77.21, 57.04, 45.85, 38.77, 33.91, 30.37, 27.66, 25.51, 23.76, 22.29, 21.03, 19.94, 18.99, 18.15, 17.41, 16.77, 16.22, 15.76};
    for (unsigned int i = 0; i<18; ++i) { mp_consts.XI0_arr[i] = XI0_ARR[i]; }


    // assign to global wrapper, ::wrapper_mp does fetch object, it isn't in global scope :: 
    wrapper_mp = mp_consts;
    
    
  } // init_wrapper_mp() functon block ends here
  
} // namespace matep block ends here

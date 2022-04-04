/*
 * This is the *.hpp file of the Strong Coupling Correction Object (SCCO)
 * 
 * Member functions are declared at here, and then defined at SCCO.cpp.
 *
 * author: Quang. Zhang (timohyva@github)
 *
 */ 

#ifndef MATEP_HPP
#define MATEP_HPP

//#include <string>
#include <iostream>
#include <cstddef>
#include <cmath>


using real_t = float;

class MATEP {
public:
        MATEP() = default; // constructor

        // interfaces of dimensional qualities
	
        real_t Tcp(real_t p); // in unit of Kelvin
        real_t Tcp_mK(real_t p); // in unit of mK
  
        real_t mEffp(real_t p);
        real_t vFp(real_t p);
        real_t xi0p(real_t p);
     // real_t N0p(real_t p);
        double N0p(real_t p);

        // interfaces of dimensionless coeficients; SC-correction parts:
  
        real_t alpha_td(real_t p, real_t T);
        real_t beta1_td(real_t p, real_t T);
        real_t beta2_td(real_t p, real_t T);
        real_t beta3_td(real_t p, real_t T);
        real_t beta4_td(real_t p, real_t T);
        real_t beta5_td(real_t p, real_t T);
        
private:
        // SI unit
        static constexpr real_t J = 1.0f, s = 1.0f, m = 1.0f, kg = 1.0f
	                        ,pi = 3.14159265358979323846264338328f;

        // physical constants for he3
        static const real_t u, m3, nm, hbar
	                    ,zeta3, c_betai ;
 
       
        // SC-data sheets arries, All associate to SCCO class
        static const real_t c1_arr[18]; 
        static const real_t c2_arr[18];
        static const real_t c3_arr[18];
        static const real_t c4_arr[18];
        static const real_t c5_arr[18];
        // ***************************************************
        // '''
        static const real_t Tc_arr[18];
        static const real_t Ms_arr[18];
        static const real_t VF_arr[18];
        static const real_t XI0_arr[18];

        // linear interpolation member:
        real_t lininterp(const real_t *cX_arr, real_t p);
};

#endif

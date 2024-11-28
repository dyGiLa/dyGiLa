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

#ifndef MATEP_HPP
#define MATEP_HPP

//#include <string>
#include <iostream>
#include <cstddef>
#include <cmath>
#include <vector>

//#include "plumbing/hila.h"

namespace matep {

 using real_t = float; // or double
  
 class Matep {
  public:

//#pragma hila loop_function
     //Matep() = default;
        //Matep(/*hila::global<matep_consts> &*/);                                      // constructor which initialize wrapper_m
        //Matep(): Switch("OFF") {};                                                  // default constructor
        //Matep(const std::string &S): Switch(S) {};                                  // metap with fudge exponent switch

  // ************************************************************************** //
  // >>>>>>>>>>>        interfaces of dimensional qualities        <<<<<<<<<<<< //
  // ************************************************************************** //

#pragma hila loop_function
        real_t Tcp(real_t p);                  // in unit of Kelvin
#pragma hila loop_function  
        real_t Tcp_mK(real_t p);               // in unit of mK
#pragma hila loop_function  
        real_t mEffp(real_t p);                // quisiparticle effective mass
#pragma hila loop_function  
        real_t vFp(real_t p);                  // Fermi velocity
#pragma hila loop_function  
        real_t xi0p(real_t p);                 // zero Temperature coherent length
#pragma hila loop_function  
        double N0p(real_t p);                  // deisty of state on Fermi surface
  

  // ************************************************************************* //
  // >>>>  interfaces of dimensionless coeficients; SC-correction parts: <<<<< //
  // ************************************************************************* //


#pragma hila loop_function
        real_t alpha_td(real_t p, real_t T);
#pragma hila loop_function  
        real_t beta1_td(real_t p, real_t T);
#pragma hila loop_function  
        real_t beta2_td(real_t p, real_t T);
#pragma hila loop_function  
        real_t beta3_td(real_t p, real_t T);
#pragma hila loop_function  
        real_t beta4_td(real_t p, real_t T);
#pragma hila loop_function  
        real_t beta5_td(real_t p, real_t T);

  // >>>>>>>>>    interfaces for beta_A, beta_B, gaps and tAB_RWS    <<<<<<<<< //

#pragma hila loop_function
        real_t beta_A_td(real_t p, real_t T);
#pragma hila loop_function  
        real_t beta_B_td(real_t p, real_t T);
#pragma hila loop_function  
        real_t gap_A_td(real_t p, real_t T);
#pragma hila loop_function  
        real_t gap_B_td(real_t p, real_t T);
#pragma hila loop_function  
        real_t gap_td(real_t p, real_t T);      // gap for given p,T, and showing message

#pragma hila loop_function  
        real_t tAB_RWS(real_t p);               

  // >>>>>>> interfaces for f_{A}, f_{B}, in unit of (1/3)(Kb Tc)^2 N(0) <<<<< //


#pragma hila loop_function
        real_t f_A_td(real_t p, real_t T);
#pragma hila loop_function  
        real_t f_B_td(real_t p, real_t T);
        

  // >>>>>>>>>>>>>>> menmber funcition Levi_Civita symbol <<<<<<<<<<<<<<<<<<<< //

#pragma hila loop_function
        real_t epsilon(int al, /*alpha*/
         	       int be, /*beta*/
		       int ga /*gamma*/);
  
 private:
  

        // linear interpolation function:
#pragma hila loop_function  
        real_t lininterp(const real_t *cX_arr, real_t p);


        // fudge expotent calculator
        // real_t exp_q(real_t p);

 }; // class Matep block ends here

} // namespace matep block ends here  

#endif

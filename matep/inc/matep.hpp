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


using real_t = float; // or double

class Matep {
public:
        Matep():
	  Switch("OFF") {};                                                  // default constructor
  /*Matep(const std::string &tempScale):
    Switch("OFF"), temperature_scale(tempScale) {};*/                    // constructor with temperature scale option  
        Matep(const std::string &S):
	  Switch(S) {};                                                      // metap with fudge exponent switch

  // ************************************************************************** //
  // >>>>>>>>>>>        interfaces of dimensional qualities        <<<<<<<<<<<< //
  // ************************************************************************** //
        real_t Fa0p(real_t p);                 //dimensionless Landau Coefficient
        real_t Tcp(real_t p);                  // in unit of Kelvin
        real_t Tcp_mK(real_t p);               // in unit of mK

        real_t tauQP(real_t p, real_t T);      // Fitted QP releaxing time in \mus
        real_t mEffp(real_t p);                // quisiparticle effective mass 
        real_t vFp(real_t p);                  // Fermi velocity
        real_t xi0p(real_t p);                 // zero Temperature coherent length
        real_t xi0GLp(real_t p);               // zero temperature GL choherent length
        real_t tGL(real_t p);                  // GL time, time unit
        double N0p(real_t p);                  // deisty of state on Fermi surface
        real_t Dd(real_t p);                   // diffusion constant/ diffusive coefficient

  // ************************************************************************* //
  // >> APIs of key properties of spheric hot bloob; SC-correction parts: <<<< //
  // ************************************************************************* //

        // time when fountier of Tc achieve maximum
        real_t t_TcMax_blob(real_t p, real_t Ttdb1, real_t Ttdb0, real_t t1);
  
        // time when fountier of Tc achieve maximum  
        real_t r_TcMax_blob(real_t p, real_t Ttdb1, real_t Ttdb0, real_t t1);

        // time when fountier of Tc shriks to vanish
        // after this T < Tc over box
        real_t t_TcVanish_blob(real_t p, real_t Ttdb1, real_t Ttdb0, real_t t1);
  

  // ************************************************************************* //
  // >>>>  interfaces of dimensionless coeficients; SC-correction parts: <<<<< //
  // ************************************************************************* //
  
        real_t alpha_td(real_t p, real_t T);
        real_t beta1_td(real_t p, real_t T);
        real_t beta2_td(real_t p, real_t T);
        real_t beta3_td(real_t p, real_t T);
        real_t beta4_td(real_t p, real_t T);
        real_t beta5_td(real_t p, real_t T);

        // damping term coefficient gamma, in unit of tGL^-1
        /* Hagen Kleinert's damping coefficient */
        real_t gamma_td(real_t p, real_t T);          

        // Quadratic H-term coefficient gz_td
        real_t gz_td(real_t p);        

  // >>>>>>>>>    interfaces for beta_A, beta_B, gaps and tAB_RWS    <<<<<<<<< //
  
        real_t beta_A_td(real_t p, real_t T);
        real_t beta_B_td(real_t p, real_t T);
        real_t gap_A_td(real_t p, real_t T);
        real_t gap_B_td(real_t p, real_t T);
        real_t gap_td(real_t p, real_t T);      // gap for given p,T, and showing message

        real_t tAB_RWS(real_t p);               

  // >>>>>>> interfaces for f_{A}, f_{B}, in unit of (1/3)(Kb Tc)^2 N(0) <<<<< //
  
        real_t f_A_td(real_t p, real_t T);
        real_t f_B_td(real_t p, real_t T);
        

  // >>>>>>>>>>>>>>> menmber funcition Levi_Civita symbol <<<<<<<<<<<<<<<<<<<< //

        real_t epsilon(int al, /*alpha*/
         	       int be, /*beta*/
		       int ga /*gamma*/);
  
private:
        // SI unit
        static constexpr real_t Kelvin = 1.0f, J = 1.0f, s = 1.0f, m = 1.0f, kg = 1.0f, mT = 1.0f
	                        ,pi = 3.14159265358979323846264338328f, E = 2.718281828459045235360287471352f
	                        ,p_pcp = 21.22f;

        // temperature switch Greywall or PLTS2000
        std::string temperature_scale;
        
        // fudge Switch, "ON" or "OFF"
        std::string Switch /*= "OFF"*/; 

        // physical constants for he3
        static const real_t u, m3, nm, mus, hbar, kb ,zeta3, c_betai, gammahbar, tau0N;
 
       
        // SC-data sheets arries, All associate to SCCO class
        static const real_t c1_arr[18]; 
        static const real_t c2_arr[18];
        static const real_t c3_arr[18];
        static const real_t c4_arr[18];
        static const real_t c5_arr[18];

        // ***************************************************
        static const real_t Fa0_arr[18];
        static const real_t Tc_arr[18];
        static const real_t Ms_arr[18];
        static const real_t VF_arr[18];
        static const real_t XI0_arr[18];

        // coefficients vector for fudge polynomial
        static const std::vector<real_t> coef4;
        // static const std::vector<real_t> coef6;

        // linear interpolation function:
        real_t lininterp(const real_t *cX_arr, real_t p);

        // fudge expotent calculator
        real_t exp_q(real_t p);
};

#endif

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

#include "matep.hpp"


//*********************************************************************
//***     member functions, interfaces of dimensional qualities     ***
//*********************************************************************

real_t
Matep::Tcp(real_t p){
  real_t Tc = lininterp(Tc_arr, p)*std::pow(10.0f,-3);
  return Tc;
}

real_t
Matep::Tcp_mK(real_t p) {return lininterp(Tc_arr, p);
}


real_t
Matep::mEffp(real_t p){  
  real_t mEff = lininterp(Ms_arr, p)*m3;;
  return mEff;
}

real_t
Matep::vFp(real_t p){
  // unit m.s^-1
  real_t vF = lininterp(VF_arr, p);
  return vF;
}

real_t
Matep::xi0p(real_t p){
  real_t xi0 = lininterp(XI0_arr, p)*nm;
  return xi0;
}  

double
Matep::N0p(real_t p){
  /*
   * the maginitude of N0p is about 10^(50), it must be double type 
   */
  double N0 = (std::pow(mEffp(p),2)*vFp(p))/((2.0f*pi*pi)*std::pow(hbar,3));
  // ((mEff(p)**(2))*vF(p))/((2*pi*pi)*(hbar**(3)))
  return N0;
}


//**********************************************************************
//***    member functions, interfaces of dimensionless coefficients  ***
//**********************************************************************

real_t
Matep::alpha_td(real_t p, real_t T){ return 1.f*(T/Tcp_mK(p)-1); }  


real_t
Matep::beta1_td(real_t p, real_t T){
  real_t beta1 = c_betai*(-1.0f + (T/Tcp_mK(p))*exp_q(p)*lininterp(c1_arr, p));

  return beta1;
}  


real_t
Matep::beta2_td(real_t p, real_t T){
  real_t beta2 = c_betai*(2.0f + (T/Tcp_mK(p))*exp_q(p)*lininterp(c2_arr, p));

  return beta2;
}  


real_t
Matep::beta3_td(real_t p, real_t T){
  real_t beta3 = c_betai*(2.0f + (T/Tcp_mK(p))*exp_q(p)*lininterp(c3_arr, p));

  return beta3;
}  


real_t
Matep::beta4_td(real_t p, real_t T){
  real_t beta4 = c_betai*(2.0f + (T/Tcp_mK(p))*exp_q(p)*lininterp(c4_arr, p));

  return beta4;
}


real_t
Matep::beta5_td(real_t p, real_t T){
  real_t beta5 = c_betai*(-2.0f + (T/Tcp_mK(p))*exp_q(p)*lininterp(c5_arr, p));

  return beta5;
}  


//**********************************************************************
//***                 beta_A, beta_B and Gaps                        ***
//**********************************************************************

real_t
Matep::beta_A_td(real_t p, real_t T){
  return beta2_td(p, T) + beta4_td(p, T) + beta5_td(p, T);
}

real_t
Matep::beta_B_td(real_t p, real_t T){
  return beta1_td(p, T) + beta2_td(p, T) + (1.f/3.f)*(beta3_td(p, T) + beta4_td(p, T) + beta5_td(p, T));
}

// A-phase gap energy, in unit of Kb * Tc
real_t
Matep::gap_A_td(real_t p, real_t T){

  if (T <= Tcp_mK(p))
    {
      real_t gap2 =-alpha_td(p, T)/(2.f*beta_A_td(p, T)); // (kb Tc)^2

      return std::sqrt(gap2);    
    }
  else //if (T > Tcp_mK(p))
    return 0.;

}

// B-phase gap energy, in unit of Kb * Tc
real_t
Matep::gap_B_td(real_t p, real_t T){

  if (T <= Tcp_mK(p))
    {
      real_t gap2 =-alpha_td(p, T)/(2.f*beta_B_td(p, T)); // (kb Tc)^2

      return std::sqrt(gap2);
    }
  else //if (T > Tcp_mK(p))
    return 0.;
}

// A general gap function with message of equlibrium phase
real_t
Matep::gap_td(real_t p, real_t T){

  if (f_A_td(p, T) > f_B_td(p, T)){
    std::cout << " \nnow p, T are: " << p << ", " << T
              << ", equlibrum bulk phase is B phase. "
              << std::endl;
    return gap_A_td(p, T);
    
  } else if (f_A_td(p, T) < f_B_td(p, T)) { 
    
    std::cout << " \nnow p, T are: " << p << ", " << T
              << ", equlibrum bulk phase is A phase. "
              << std::endl;
    return gap_B_td(p, T);   

  } else {

    if (
	// (f_A_td(p, T) == f_B_td(p, T))
	// && (T < Tcp_mK(p))
	T < Tcp_mK(p)
       ){

       std::cout << " \nnow p, T are: " << p << ", " << T
                 << ", and A and B degenerate, return as -1. "
                 << std::endl;
       return -1.f;

    } else  //(
        	// (f_A_td(p, T) == f_B_td(p, T))
	        // && (T > Tcp_mK(p))
            // )
	  {

            std::cout << " \nnow p, T are: " << p << ", " << T
	            << ", system is in normal phase. "
	            << std::endl;
	    return 0.f;

    }

  }
}

// tAB_RWS 2019
real_t
Matep::tAB_RWS(real_t p){
  real_t t = 1.f/(3.f*lininterp(c1_arr, p)
           + lininterp(c3_arr, p)
           - 2.f*lininterp(c4_arr, p)
           - 2.f*lininterp(c5_arr, p));
  return t;
}


// A-Phase free energy density in unit of (1/3)(Kb Tc)^2 N(0)
real_t
Matep::f_A_td(real_t p, real_t T)
{
  if (T <= Tcp_mK(p))
    {
     return (-1.f/4.f)*(std::pow(alpha_td(p, T),2.f))/beta_A_td(p, T);    
    }
  else //if (T > Tcp_mK(p))
    return 0.;
    
}

// B-Phase free energy density in unit of (1/3)(Kb Tc)^2 N(0)
real_t
Matep::f_B_td(real_t p, real_t T)
{
  if (T <= Tcp_mK(p))
    {
     return (-1.f/4.f)*(std::pow(alpha_td(p, T),2.f))/beta_B_td(p, T);
    }
  else //if (T > Tcp_mK(p))
    return 0.;

}


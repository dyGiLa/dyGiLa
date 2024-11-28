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
#include "matep.hpp"

//#include "plumbing/hila.h"

//#include "matep_consts.hpp"
#include "plumbing/globals.h"

// declearation of hila::global<matep::matep_consts> wrapper_mp
// the initialization is done by matep::init_wrapper_mp() call,
// decleared in matep_namespace_utils.hpp
hila::global<matep::matep_consts> wrapper_mp; 


//*********************************************************************
//***     member functions, interfaces of dimensional qualities     ***
//*********************************************************************

namespace matep {

real_t
Matep::Tcp(real_t p){
  real_t Tc = lininterp(wrapper_mp().Tc_arr, p)*(10.0e-3f);
  //real_t Tc = p*p + p*p*p;
  return Tc;
}

real_t
Matep::Tcp_mK(real_t p) {return lininterp(wrapper_mp().Tc_arr, p);
}


real_t
Matep::mEffp(real_t p){  
  real_t mEff = lininterp(wrapper_mp().Ms_arr, p)*wrapper_mp().m3;
  //real_t mEff = p*p + p*p*p;  
  return mEff;
}

real_t
Matep::vFp(real_t p){
  // unit m.s^-1
  real_t vF = lininterp(wrapper_mp().VF_arr, p);
  //real_t vF = p*p + p*p*p;  
  return vF;
}

real_t
Matep::xi0p(real_t p){
  real_t xi0 = lininterp(wrapper_mp().XI0_arr, p)*wrapper_mp().nm;
  //real_t xi0 = p*p + p*p*p;  
  return xi0;
}  

double
Matep::N0p(real_t p){
  /*
   * the maginitude of N0p is about 10^(50), it must be double type 
   */
  double N0 = (mEffp(p)*mEffp(p)*vFp(p))/((2.0f*wrapper_mp().pi*wrapper_mp().pi)*(wrapper_mp().hbar*wrapper_mp().hbar*wrapper_mp().hbar));
  //real_t N0 = p*p + p*p*p;  
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
  real_t beta1 = wrapper_mp().c_betai*(-1.0f + (T/Tcp_mK(p))*lininterp(wrapper_mp().c1_arr, p));
  return beta1;
}  


real_t
Matep::beta2_td(real_t p, real_t T){
  real_t beta2 = wrapper_mp().c_betai*(2.0f + (T/Tcp_mK(p))*lininterp(wrapper_mp().c2_arr, p));

  return beta2;
}  


real_t
Matep::beta3_td(real_t p, real_t T){
  real_t beta3 = wrapper_mp().c_betai*(2.0f + (T/Tcp_mK(p))*lininterp(wrapper_mp().c3_arr, p));

  return beta3;
}  


real_t
Matep::beta4_td(real_t p, real_t T){
  real_t beta4 = wrapper_mp().c_betai*(2.0f + (T/Tcp_mK(p))*lininterp(wrapper_mp().c4_arr, p));

  return beta4;
}


real_t
Matep::beta5_td(real_t p, real_t T){
  real_t beta5 = wrapper_mp().c_betai*(-2.0f + (T/Tcp_mK(p))*lininterp(wrapper_mp().c5_arr, p));

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
    // hila::out << " \nnow p, T are: " << p << ", " << T
    //           << ", equlibrum bulk phase is B phase. "
    //           << std::endl;
    return gap_A_td(p, T);
    
  } else if (f_A_td(p, T) < f_B_td(p, T)) { 
    
    // hila::out << " \nnow p, T are: " << p << ", " << T
    //           << ", equlibrum bulk phase is A phase. "
    //           << std::endl;
    return gap_B_td(p, T);   

  } else {

    if (
	// (f_A_td(p, T) == f_B_td(p, T))
	// && (T < Tcp_mK(p))
	T < Tcp_mK(p)
       ){

       // hila::out << " \nnow p, T are: " << p << ", " << T
       //           << ", and A and B degenerate, return as -1. "
       //           << std::endl;
       return -1.f;

    } else  //(
        	// (f_A_td(p, T) == f_B_td(p, T))
	        // && (T > Tcp_mK(p))
            // )
	  {

            // hila::out << " \nnow p, T are: " << p << ", " << T
	    //         << ", system is in normal phase. "
	    //         << std::endl;
	    return 0.f;

    }

  }
}

// tAB_RWS 2019
real_t
Matep::tAB_RWS(real_t p){
  real_t t = 1.f/(3.f*lininterp(wrapper_mp().c1_arr, p)
           + lininterp(wrapper_mp().c3_arr, p)
           - 2.f*lininterp(wrapper_mp().c4_arr, p)
           - 2.f*lininterp(wrapper_mp().c5_arr, p));
  // real_t t = p*p*p + p*p;
  return t;
}


// A-Phase free energy density in unit of (1/3)(Kb Tc)^2 N(0)
real_t
Matep::f_A_td(real_t p, real_t T)
{
  if (T <= Tcp_mK(p))
    {
     return (-1.f/4.f)*((alpha_td(p, T)*alpha_td(p, T)))/beta_A_td(p, T);    
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
     return (-1.f/4.f)*((alpha_td(p, T)*alpha_td(p, T)))/beta_B_td(p, T);
    }
  else //if (T > Tcp_mK(p))
    return 0.;

}

} // namespace matep block ends here

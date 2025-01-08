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
Matep::Fa0p(real_t p) {return lininterp(Fa0_arr, p);
}

real_t
Matep::Tcp(real_t p){
  real_t Tc = lininterp(Tc_arr, p)*std::pow(10.0f,-3);
  return Tc;
}

real_t
Matep::Tcp_mK(real_t p) {return lininterp(Tc_arr, p);
}

real_t
Matep::tauQP(real_t p, real_t T){

  real_t t = T/Tcp_mK(p);
  
  return (-2.8277653941154503e6
	  + 6.284350363777282e7*t - 6.286343309406481e8*std::pow(t,2) + 3.7463162643827033e9*std::pow(t,3)
	  - 1.4784434682274881e10*std::pow(t,4) + 4.042411599324316e10*std::pow(t,5) - 7.713192342581464e10*std::pow(t,6)
	  + 9.84181184085264e10*std::pow(t,7) - 6.836431172791403e10*std::pow(t,8) - 1.4180711852241095e10*std::pow(t,9)
	  + 9.61270986677446e10*std::pow(t,10) - 1.2002060849381773e11*std::pow(t,11) + 8.539279623542038e10*std::pow(t,12)
	  - 3.750912433469162e10*std::pow(t,13) + 9.530297512559393e9*std::pow(t,14) - 1.0790099723106315e9*std::pow(t,15));

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

real_t
Matep::xi0GLp(real_t p){ return xi0p(p)*std::sqrt((7.f*zeta3)/20.f); }

real_t
Matep::tGL(real_t p){ return 1.290994449*(xi0GLp(p)/vFp(p)); }

double
Matep::N0p(real_t p){
  /*
   * the maginitude of N0p is about 10^(50), it must be double type 
   */
  double N0 = (std::pow(mEffp(p),2)*vFp(p))/((2.0f*pi*pi)*std::pow(hbar,3));
  // ((mEff(p)**(2))*vF(p))/((2*pi*pi)*(hbar**(3)))
  return N0;
}

real_t
Matep::Dd(real_t p){
  real_t vF = vFp(p);
  real_t xi0GL = xi0GLp(p);

  // convert ratio is made from SI unit value number
  real_t tGLxiGL2_ratio = tGL(p)/(xi0GL * xi0GL);
  // return diffussion constant in unit of xiGL^2.tGL^-1
  return vF * vF * tau0N * tGLxiGL2_ratio;
}

// ************************************************************************* //
// >> APIs of key properties of spheric hot bloob; SC-correction parts: <<<< //
// ************************************************************************* //


real_t
Matep::t_TcMax_blob(real_t p, real_t Ttdb1, real_t Ttdb0, real_t t1) {

  // t1 could be in unit of tGL, Tx has same 
  real_t TcpmK  = Tcp_mK(p);
  real_t Tx     = (Ttdb1 - Ttdb0) * TcpmK;
  real_t T0     = Ttdb0 * TcpmK;
  
  return (std::pow(-1, 0.6666666666666666) * t1 * std::pow(Tx, 0.6666666666666666))
         /(E * std::pow((T0 - TcpmK), 0.6666666666666666));
}

real_t
Matep::r_TcMax_blob(real_t p, real_t Ttdb1, real_t Ttdb0, real_t t1) {

  real_t TcpmK  = Tcp_mK(p);  
  real_t Tx     = (Ttdb1 - Ttdb0) * TcpmK;
  real_t T0     = Ttdb0 * TcpmK;
  
  return sqrt(6./E) * sqrt(Dd(p) * t1)* std::pow((Tx/(-T0 + TcpmK)),0.3333333333333333);
}

real_t
Matep::t_TcVanish_blob(real_t p, real_t Ttdb1, real_t Ttdb0, real_t t1) {
  real_t TcpmK = Tcp_mK(p);  
  real_t Tx    = (Ttdb1 - Ttdb0) * TcpmK;
  real_t T0    = Ttdb0 * TcpmK;
  
  return (std::pow(-1., 0.6666666666666666) * t1 * std::pow(Tx, 0.6666666666666666))
         /std::pow((T0 - TcpmK), 0.6666666666666666);
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

real_t
Matep::gz_td(real_t p){
  real_t gz = 5.f*c_betai
              *(1./((1+Fa0p(p))*(1+Fa0p(p))))
              *(gammahbar/(kb*Tcp(p)))
              *(gammahbar/(kb*Tcp(p)));

  return gz;
}

// real_t
// Matep::gamma_td(real_t p, real_t T){ return tGL(p)/(tauQP(p, T)*mus); }  

real_t
Matep::gamma_td(real_t p, real_t T, real_t pM){
  real_t gamma_C = ((pi*pi)/4) * sqrt(3./5.) * sqrt(20./(7.*zeta3))
         ,gtd;
  real_t t = T/Tcp_mK(p);
  
  if (t >= 1.0 /*T >= Tcp_mK(p)*/)
     gtd = 1.;
  else
    {
      if ( pM == 4.0 )
	gtd = std::pow(t, 4.); // A-Phase
      else if ( pM == 2.0f )
	gtd = exp(3.349621654292973*(1 - 1./t)); // B-phase
      else if ( pM == 1.0f )
	gtd = std::pow(t, 4.); // planar-phase
      else if ( pM == 3.0f )
	gtd = std::pow(t, 4.); // polar-phase
      else if ( pM == 0.0f )
	gtd = std::pow(t, 4.);		
    }

  return gtd * gamma_C;
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


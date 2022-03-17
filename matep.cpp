/*
 * This is the *.cpp file of the Strong Coupling Correction Object (SCCO)
 * 
 * Member functions are declared at SCCO.hpp, and then defined at here.
 *
 * author: Quang. Zhang (timohyva@github)
 *
 */ 


#include <iostream>
#include <cstddef>
#include <cmath>

#include "matep.hpp"

//********************************************************************
//********************************************************************
//''' data sheets of  strong coupling corrected material parameters
//'

//constexpr long double SCCO::pi = = 3.14159265358979323846264338328L;
const real_t MATEP::c1_arr[18] = {-0.0098, -0.0127, -0.0155, -0.0181, -0.0207, -0.0231, -0.0254, -0.0275, -0.0295, -0.0314, -0.0330, -0.0345, -0.0358, -0.0370, -0.0381, -0.0391, -0.0402, -0.0413};
const real_t MATEP::c2_arr[18] = {-0.0419, -0.0490, -0.0562, -0.0636, -0.0711, -0.0786, -0.0861, -0.0936, -0.1011, -0.1086, -0.1160, -0.1233, -0.1306, -0.1378, -0.1448, -0.1517, -0.1583, -0.1645};
const real_t MATEP::c3_arr[18] = {-0.0132, -0.0161, -0.0184, -0.0202, -0.0216, -0.0226, -0.0233, -0.0239, -0.0243, -0.0247, -0.0249, -0.0252, -0.0255, -0.0258, -0.0262, -0.0265, -0.0267, -0.0268};
const real_t MATEP::c4_arr[18] = {-0.0047, -0.0276, -0.0514, -0.0760, -0.1010, -0.1260, -0.1508, -0.1751, -0.1985, -0.2208, -0.2419, -0.2614, -0.2795, -0.2961, -0.3114, -0.3255, -0.3388, -0.3518};
const real_t MATEP::c5_arr[18] = {-0.0899, -0.1277, -0.1602, -0.1880, -0.2119, -0.2324, -0.2503, -0.2660, -0.2801, -0.2930, -0.3051, -0.3167, -0.3280, -0.3392, -0.3502, -0.3611, -0.3717, -0.3815};

const real_t MATEP::Tc_arr[18] = {0.929, 1.181, 1.388, 1.560, 1.705, 1.828, 1.934, 2.026, 2.106, 2.177, 2.239, 2.293, 2.339, 2.378, 2.411, 2.438, 2.463, 2.486}; // mK
const real_t MATEP::Ms_arr[18] = {2.80, 3.05, 3.27, 3.48, 3.68, 3.86, 4.03, 4.20, 4.37, 4.53, 4.70, 4.86, 5.02, 5.18, 5.34, 5.50, 5.66, 5.82}; // in unit of helium-3 atom
const real_t MATEP::VF_arr[18] = {59.03, 55.41, 52.36, 49.77, 47.56, 45.66, 44.00, 42.51, 41.17, 39.92, 38.74, 37.61, 36.53, 35.50, 34.53, 33.63, 32.85, 32.23}; // fermi velosity, m.s^-1
const real_t MATEP::XI0_arr[18] = {77.21, 57.04, 45.85, 38.77, 33.91, 30.37, 27.66, 25.51, 23.76, 22.29, 21.03, 19.94, 18.99, 18.15, 17.41, 16.77, 16.22, 15.76};


//*********************************************************************
//*********************************************************************
//'''     member functions, interfaces of dimensional qualities
//'

real_t
MATEP::Tcp(real_t p){
  real_t Tc = lininterp(Tc_arr, p)*std::pow(10.0f,-3);
  return Tc;
}

real_t
MATEP::mEffp(real_t p){
  constexpr float kg = 1.0f;
  const float u = 1.66053906660f*(std::pow(10.0f,-27))*kg;
  const float m3 = 3.016293f*u;
  
  real_t mEff = lininterp(Ms_arr, p)*m3;;
  return mEff;
}

real_t
MATEP::vFp(real_t p){
  // unit m.s^-1
  real_t vF = lininterp(VF_arr, p);
  return vF;
}

real_t
MATEP::xi0p(real_t p){
  constexpr real_t m = 1.0f;
  const real_t nm = (std::pow(10.0f, -9))*m; 
  real_t xi0 = lininterp(XI0_arr, p)*nm;
  return xi0;
}  

real_t
MATEP::N0p(real_t p){
  constexpr float J = 1.0f, s = 1.0f, pi = 3.14159265358979323846264338328f;;
  const float hbar = 1.054571817f*(std::pow(10.0f,-34))*J*s;
  real_t N0 = (std::pow(mEffp(p),2)*vFp(p))/((2.0f*pi*pi)*std::pow(hbar,3));
  // ((mEff(p)**(2))*vF(p))/((2*pi*pi)*(hbar**(3)))
  return N0;
}


//**********************************************************************
//**********************************************************************
//'''     member functions, interfaces of dimensionless coefficients

real_t
MATEP::alpha_bar(real_t p, real_t T){ return (1.f/3.f)*(T/Tcp(p)-1); }  


real_t
MATEP::beta1_bar(real_t p, real_t T){
  constexpr float pi = 3.14159265358979323846264338328f;
  const float zeta3 = std::riemann_zeta(3.0f);
  const float c_betai = (7.0f*zeta3)/(240.0f*pi*pi);
  real_t beta1 = c_betai*(-1.0f + (T/Tcp(p))*lininterp(c1_arr, p));

  return beta1;
}  


real_t
MATEP::beta2_bar(real_t p, real_t T){
  constexpr float pi = 3.14159265358979323846264338328f;
  const float zeta3 = std::riemann_zeta(3.0f);
  const float c_betai = (7.0L*zeta3)/(240.0L*pi*pi);
  real_t beta1 = c_betai*(2.0f + (T/Tcp(p))*lininterp(c2_arr, p));

  return beta1;
}  


real_t
MATEP::beta3_bar(real_t p, real_t T){
  constexpr float pi = 3.14159265358979323846264338328f;
  const float zeta3 = std::riemann_zeta(3.0f);
  const float c_betai = (7.0L*zeta3)/(240.0L*pi*pi);
  real_t beta1 = c_betai*(2.0f + (T/Tcp(p))*lininterp(c3_arr, p));

  return beta1;
}  


real_t
MATEP::beta4_bar(real_t p, real_t T){
  constexpr float pi = 3.14159265358979323846264338328f;
  const float zeta3 = std::riemann_zeta(3.0f);
  const float c_betai = (7.0f*zeta3)/(240.0f*pi*pi);
  real_t beta1 = c_betai*(2.0f + (T/Tcp(p))*lininterp(c4_arr, p));

  return beta1;
}


real_t
MATEP::beta5_bar(real_t p, real_t T){
  constexpr float pi = 3.14159265358979323846264338328f;
  const float zeta3 = std::riemann_zeta(3.0f);
  const float c_betai = (7.0f*zeta3)/(240.0f*pi*pi);
  real_t beta1 = c_betai*(-2.0f + (T/Tcp(p))*lininterp(c5_arr, p));

  return beta1;
}  



//**********************************************************************
//**********************************************************************
//'''                  linear intepolation function
//'

real_t
MATEP::lininterp(const real_t *cX_arr, real_t p){
  float pk, pk1, fp;
  size_t k, k1;

  if ((p >= 0.0) && (p < 2.0)) { pk = 0.0f; k = 0; }

  if ((p >= 2.0) && (p < 4.0)) { pk = 2.0f; k = 1; }
 
  if ((p >= 4.0) && (p < 6.0)) { pk = 4.0f; k = 2; }

  if ((p >= 6.0) && (p < 8.0)) { pk = 6.0f; k = 3; }

  if ((p >= 8.0) && (p < 10.0)) { pk = 8.0f; k = 4; }

  if ((p >= 10.0) && (p < 12.0)) { pk = 10.0f; k = 5; }

  if ((p >= 12.0) && (p < 14.0)) { pk = 12.0f; k = 6; }

  if ((p >= 14.0) && (p < 16.0)) { pk = 14.0f; k = 7; }

  if ((p >= 16.0) && (p < 18.0)) { pk = 16.0f; k = 8; }

  if ((p >= 18.0) && (p < 20.0)) { pk = 18.0f; k = 9; }

  if ((p >= 20.0) && (p < 22.0)) { pk = 20.0f; k = 10; }

  if ((p >= 22.0) && (p < 24.0)) { pk = 22.0f; k = 11; }

  if ((p >= 24.0) && (p < 26.0)) { pk = 24.0f; k = 12; }

  if ((p >= 26.0) && (p < 28.0)) { pk = 26.0f; k = 13; }

  if ((p >= 28.0) && (p < 30.0)) { pk = 28.0f; k = 14; }

  if ((p >= 30.0) && (p < 32.0)) { pk = 30.0f; k = 15; }

  if ((p >= 32.0) && (p < 34.0)) { pk = 32.0f; k = 16; }

  fp = ((cX_arr[k+1]-cX_arr[k])/2.0)*(p-pk)+cX_arr[k];
  return fp; 
}










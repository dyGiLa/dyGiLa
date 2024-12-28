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


//********************************************************************
//***          static physical constants of he3 members            ***
//********************************************************************
const real_t Matep::kb = 1.380649*(std::pow(10.0f, -23))*J*1.0f;
const real_t Matep::u = 1.66053906660f*(std::pow(10.0f,-27))*kg;
const real_t Matep::m3 = 3.016293f*u;
const real_t Matep::nm = 1.0e-9F; // (std::pow(10.0f, -9))*m;
const real_t Matep::mus = 1e-6*s;
const real_t Matep::hbar = 1.054571817f*(std::pow(10.0f,-34))*J*s;
const real_t Matep::zeta3 = std::riemann_zeta(3.0f);
const real_t Matep::c_betai = (7.0f*zeta3)/(80.0f*pi*pi);
const real_t Matep::gammahbar = -34.2040866*(std::pow(10.0f,-31))*J*1.0f; // in unit of J-mT^-1
const real_t Matep::tau0N = 0.3e-6*s; // normal phase QP life time 0.3 mus in unit of s


// *******************************************************************
// >> data sheets of  strong coupling corrected material parameters <<
//'*******************************************************************

const real_t Matep::c1_arr[18] = {-0.0098, -0.0127, -0.0155, -0.0181, -0.0207, -0.0231, -0.0254, -0.0275, -0.0295, -0.0314, -0.0330, -0.0345, -0.0358, -0.0370, -0.0381, -0.0391, -0.0402, -0.0413};
const real_t Matep::c2_arr[18] = {-0.0419, -0.0490, -0.0562, -0.0636, -0.0711, -0.0786, -0.0861, -0.0936, -0.1011, -0.1086, -0.1160, -0.1233, -0.1306, -0.1378, -0.1448, -0.1517, -0.1583, -0.1645};
const real_t Matep::c3_arr[18] = {-0.0132, -0.0161, -0.0184, -0.0202, -0.0216, -0.0226, -0.0233, -0.0239, -0.0243, -0.0247, -0.0249, -0.0252, -0.0255, -0.0258, -0.0262, -0.0265, -0.0267, -0.0268};
const real_t Matep::c4_arr[18] = {-0.0047, -0.0276, -0.0514, -0.0760, -0.1010, -0.1260, -0.1508, -0.1751, -0.1985, -0.2208, -0.2419, -0.2614, -0.2795, -0.2961, -0.3114, -0.3255, -0.3388, -0.3518};
const real_t Matep::c5_arr[18] = {-0.0899, -0.1277, -0.1602, -0.1880, -0.2119, -0.2324, -0.2503, -0.2660, -0.2801, -0.2930, -0.3051, -0.3167, -0.3280, -0.3392, -0.3502, -0.3611, -0.3717, -0.3815};

const real_t Matep::Fa0_arr[18] = {-0.7226, -0.7317, -0.7392, -0.7453, -0.7503, -0.7544, -0.7580, -0.7610, -0.7637, -0.7661, -0.7684, -0.7705, -0.7725, -0.7743, -0.7758, -0.7769, -0.7775, -0.7775}; 
const real_t Matep::Tc_arr[18] = {0.929, 1.181, 1.388, 1.560, 1.705, 1.828, 1.934, 2.026, 2.106, 2.177, 2.239, 2.293, 2.339, 2.378, 2.411, 2.438, 2.463, 2.486}; // mK
const real_t Matep::Ms_arr[18] = {2.80, 3.05, 3.27, 3.48, 3.68, 3.86, 4.03, 4.20, 4.37, 4.53, 4.70, 4.86, 5.02, 5.18, 5.34, 5.50, 5.66, 5.82}; // in unit of helium-3 atom
const real_t Matep::VF_arr[18] = {59.03, 55.41, 52.36, 49.77, 47.56, 45.66, 44.00, 42.51, 41.17, 39.92, 38.74, 37.61, 36.53, 35.50, 34.53, 33.63, 32.85, 32.23}; // fermi velosity, m.s^-1
const real_t Matep::XI0_arr[18] = {77.21, 57.04, 45.85, 38.77, 33.91, 30.37, 27.66, 25.51, 23.76, 22.29, 21.03, 19.94, 18.99, 18.15, 17.41, 16.77, 16.22, 15.76}; // nm 1e-9 m


//********************************************************************
// ***          data sheet of fudge exponet polynomial             ***
//********************************************************************


const std::vector<real_t> Matep::coef4 = {-6.00498973e-03, -1.01758101e-02,  1.46969023e-03, -1.14870022e-04, 4.11400719e-06};


//const std::vector<real_t> Matep::coef4 = {-6.00498973*std::pow(10.f,-3), -1.01758101*std::pow(10.f,-2), 1.46969023*std::pow(10.f,-3), -1.14870022*std::pow(10.f,-4), 4.11400719*std::pow(10.f,-6)};


//**********************************************************************
//***              public member: Levi-Civita symbol                 ***
//**********************************************************************

real_t
Matep::epsilon(int al, int be, int ga)
{
  if (
      (al == 0 && be == 1 && ga == 2)   // 123
      ||(al == 1 && be == 2 && ga == 0) // 231
      ||(al == 2 && be == 0 && ga == 1) // 312
     )
    {return 1.0;}
  else if
     (
      (al == 2 && be == 1 && ga == 0)   // 321
      ||(al == 0 && be == 2 && ga == 1) // 132
      ||(al == 1 && be == 0 && ga == 2) // 213
     )
    {return -1.0;}
  else if ((al == be) || (al == ga) || (be == ga))
    {return 0.0;}
  else
    {return 0.0;}
}  


//**********************************************************************
//***              private method : fudge expotential                ***
//**********************************************************************

real_t
Matep::exp_q(real_t p){
  // 4th-order polynomial of q for Greywall scale

  real_t q = 0.f, defp_G;
  defp_G = p - p_pcp;
  
  // std::cout << " \n defp_G is " << defp_G << std::endl;
  if (Switch == "ON") {
    
    if (defp_G >= 0.f){
      
      for (unsigned co = 0; co < 5; ++co){
      	q += coef4[co]*(std::pow(defp_G,co));
      }
    } else {
      q = 0.f;
    } 
    
    return std::exp(q);

  } else if (Switch == "OFF") {
    q = 0.f;
    //std::cout << " q is " << q << std::endl;
    return std::exp(0.f);

  } else {

    std::cout << " \n Switch must be \"ON\" or \"OFF\"! " << std::endl;
    return 1;
  }
  
}  

//**********************************************************************
//***       private method :  linear intepolation function           ***
//**********************************************************************

real_t
Matep::lininterp(const real_t *cX_arr, real_t p){
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


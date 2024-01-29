#include <iostream>
//#include <cstddef>
//#include <cmath>

#include "matep.hpp"

using real_t = float;
using counter = int;

int main(){

  real_t kb = 1.380649*(std::pow(10.0f, -23));

  // real_t T[6] = {2.3*std::pow(10.f,-3),0.2f*std::pow(10.f,-3),1.21f*std::pow(10.f,-3),0.456f*std::pow(10.f,-3),1.5f*std::pow(10.f,-3),1.9f*std::pow(10.f,-3)}; //Kelvin
  real_t T[9] = {0.5f, 2.f, 2.3f, 0.2f, 1.21, 0.456f, 1.5f, 1.9f, 3.0f}; // mK
  real_t p[9] = {15.f, 20.f, 26.f, 32.f, 21.50f, 23.9f, 27.5f, 33.f, 20.f}; //bar
  
  Matep MP; // object with with fudge switch "OFF"
  //Matep MP("ON"); // object with fudge switch "ON", constructor overload std::string Switch

  

  for (counter i = 0; i <= 8; ++i){

    std::cout << "\n***************************************************\n "
              << "\n***************************************************\n "
	      << std::endl;
    
    // std::cout << "\n pressure is: " << p[i] << "bar" << " and Tc is " << MP.Tcp(p[i]) << "K" << std::endl;
    std::cout << "\n pressure is: " << p[i] << "bar" <<", temperature is " << T[i] << "mK" << " and Tc is " << MP.Tcp_mK(p[i]) << "mK" << std::endl;
    std::cout << "\n density of state is " << MP.N0p(p[i]) << "J^(-1).m^(-3), and effective mass of quisiparticle is " << MP.mEffp(p[i]) << "kg" <<std::endl;
    std::cout << "\n 0 tempreture coherent length is " << MP.xi0p(p[i]) << "m^(-3), and Fermi velocity is " <<MP.vFp(p[i]) << "m.s^(-1)" << std::endl;

    std::cout << " \n " << std::endl;

    std::cout << "beta1_td = " << MP.beta1_td(p[i], T[i]) << std::endl;
    std::cout << "beta2_td = " << MP.beta2_td(p[i], T[i]) << std::endl;
    std::cout << "beta3_td = " << MP.beta3_td(p[i], T[i]) << std::endl;
    std::cout << "beta4_td = " << MP.beta4_td(p[i], T[i]) << std::endl;
    std::cout << "beta5_td = " << MP.beta5_td(p[i], T[i]) << std::endl;
    std::cout << "alpha_td = " << MP.alpha_td(p[i], T[i]) << "\n\n===================== " << std::endl;

    // std::cout << "\n gap_A = " << MP.gap_A_td(p[i], T[i])*(MP.Tcp(p[i]))*kb << std::endl;
    // std::cout << "\n gap_B = " << MP.gap_B_td(p[i], T[i])*(MP.Tcp(p[i]))*kb << std::endl;

    std::cout << "\n gap_A = " << MP.gap_A_td(p[i], T[i]) << std::endl;
    std::cout << "\n gap_B = " << MP.gap_B_td(p[i], T[i]) << std::endl;
    std::cout << "\n gap_td returns " << MP.gap_td(p[i], T[i]) << std::endl;

    std::cout << "\n f_A = " << MP.f_A_td(p[i], T[i]) << std::endl;
    std::cout << "\n f_B = " << MP.f_B_td(p[i], T[i]) << std::endl;

    
    //MP.read_coef();

    std::cout <<  "\n===================== " << std::endl;

    real_t gap = MP.gap_td(p[i], T[i]);
    std::cout << "And the gap is: " << gap << " kb*Tc(p) "<< std::endl;

    std::cout << "\n***************************************************\n "
              << "\n***************************************************\n "
	      << std::endl;
    std::cout << "\n >>>--------------------------------------------------------------------------------<<<\n "
              << std::endl;

  }

  return 0;


}  

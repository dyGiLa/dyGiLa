#include <iostream>
//#include <cstddef>
//#include <cmath>

#include "matep.hpp"

using real_t = float;
using counter = int;

int main(){

  real_t kb = 1.380649*(std::pow(10.0f, -23));

  // real_t T[6] = {2.3*std::pow(10.f,-3),0.2f*std::pow(10.f,-3),1.21f*std::pow(10.f,-3),0.456f*std::pow(10.f,-3),1.5f*std::pow(10.f,-3),1.9f*std::pow(10.f,-3)}; //Kelvin
  real_t T[6] = {2.3f, 0.2f, 1.21, 0.456f, 1.5f, 1.9f}; // mK
  real_t p[6] = {26.f, 32.f, 21.50f, 23.9f, 27.5f, 33.f}; //bar
  
  MATEP MP; // object with with fudge switch "OFF"
  //MATEP MP("ON"); // object with fudge switch "ON", constructor overload std::string Switch

  

  for (counter i = 0; i <= 5; ++i){

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


    //MP.read_coef();

    std::cout <<  "\n===================== " << " \n\n " << std::endl;
  }

  return 0;


}  

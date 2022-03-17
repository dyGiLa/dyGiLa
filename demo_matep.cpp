#include <iostream>
//#include <cstddef>
//#include <cmath>

#include "matep.hpp"

using real_t = float;
using counter = int;

int main(){

  real_t T[5] = {0.2f*std::pow(10.f,-3),1.21f*std::pow(10.f,-3),0.456f*std::pow(10.f,-3),1.5f*std::pow(10.f,-3),1.9f*std::pow(10.f,-3)}; //Kelvin
  real_t p[5] = {32.f, 17.f, 0.11f, 26.f, 33.f}; //bar
  
  MATEP MP;

  

  for (counter i = 0; i < 5; ++i){

    std::cout << "pressure is: " << p[i] << "bar" << " and Tc is " << MP.Tcp(p[i]) << "K" << std::endl;
    // std::cout << something << " " <<  MP.alpha_bar(p, T) << std::endl;
    // std::cout << MP.N0p(p) << std::endl;

    std::cout << " \n " << std::endl;

    std::cout << "beta1_bar = " << MP.beta1_bar(p[i], T[i]) << std::endl;
    std::cout << "beta2_bar = " << MP.beta2_bar(p[i], T[i]) << std::endl;
    std::cout << "beta3_bar = " << MP.beta3_bar(p[i], T[i]) << std::endl;
    std::cout << "beta4_bar = " << MP.beta4_bar(p[i], T[i]) << std::endl;
    std::cout << "beta5_bar = " << MP.beta5_bar(p[i], T[i]) << std::endl;
    std::cout << "alpha_bar = " << MP.alpha_bar(p[i], T[i]) << std::endl;
  }

  return 0;


}  

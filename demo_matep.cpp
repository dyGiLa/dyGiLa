
#include <iostream>
//#include <cstddef>
//#include <cmath>

#include "matep.hpp"


int main(){

  long double T = 1.5*std::pow(10,-3);
  long double p = 32.0;
  long double something = 3.14159265358979323846264338328L;

  MATEP MP;

  std::cout << MP.beta1_bar(p, T) << " and Tc is " << MP.Tcp(p) << std::endl;
  std::cout << something << MP.alpha_bar(p, T) << std::endl;
  std::cout << MP.N0p(p) << std::endl;

  return 0;


}  

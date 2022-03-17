### He3-simulator
Code for simulating the order parameter of superfluid Helium 3.

Cloned from onsim (Asier Lopez-Eiguren).

Uses HILA framework (Kari Rummukainen et al): https://bitbucket.org/Kari_Rummukainen/hila/src/master/

Create branch materialp for material parameter (Quang Zhang, timohyva@github)

### demo of matep with five (p, T) pairs
~~~ C++
#include <iostream>
#include "matep.hpp"

using real_t = float;
using counter = int;

int main(){

  real_t T[5] = {0.2f*std::pow(10.f,-3),1.21f*std::pow(10.f,-3),0.456f*std::pow(10.f,-3),1.5f*std::pow(10.f,-3),1.9f*std::pow(10.f,-3)}; //Kelvin
  real_t p[5] = {32.f, 17.f, 0.11f, 26.f, 33.f}; //bar
  
  MATEP MP;
  for (counter i = 0; i < 5; ++i){

    std::cout << "pressure is: " << p[i] << "bar" << " and Tc is " << MP.Tcp(p[i]) << "K" << std::endl;
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

~~~
compiling & run the demo_matep.cpp
~~~ bash
 g++ -std=c++2a -g -o demo.app demo_matep.cpp matep.cpp matep.hpp && ./demo.app
~~~
out put
~~~bash

pressure is: 32bar and Tc is 0.002463K
 
 
beta1_bar = -0.00356392
beta2_bar = 0.00705898
beta3_bar = 0.00709694
beta4_bar = 0.00700691
beta5_bar = -0.00721186
alpha_bar = -0.306266
...
...
...

~~~
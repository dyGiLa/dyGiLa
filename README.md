### He3-simulator
Code for simulating the order parameter of superfluid Helium 3.

Cloned from onsim (Asier Lopez-Eiguren).

Uses HILA framework (Kari Rummukainen et al): https://bitbucket.org/Kari_Rummukainen/hila/src/master/

Create branch materialp for material parameter (Quang Zhang, timohyva@github)

### demo of matep with five (p, T) pairs
~~~ C++
#include <iostream>
//#include <cstddef>
//#include <cmath>

#include "matep.hpp"

using real_t = float;
using counter = int;

int main(){

  // real_t T[6] = {2.3*std::pow(10.f,-3),0.2f*std::pow(10.f,-3),1.21f*std::pow(10.f,-3),0.456f*std::pow(10.f,-3),1.5f*std::pow(10.f,-3),1.9f*std::pow(10.f,-3)}; //Kelvin
  real_t T[6] = {2.3f, 0.2f, 1.21, 0.456f, 1.5f, 1.9f}; // mK
  real_t p[6] = {26.f, 32.f, 17.f, 0.11f, 26.f, 33.f}; //bar
  
  MATEP MP;

  for (counter i = 0; i < 5; ++i){

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

 pressure is: 17bar, temperature is 1.21mK and Tc is 2.1415mK

 density of state is 8.7005e+50J^(-1).m^(-3), and effective mass of quisiparticle is 2.22886e-26kg

 0 tempreture coherent length is 2.3025e-08m^(-3), and Fermi velocity is 40.545m.s^(-1)
 
 
beta1_td = -0.0108403
beta2_td = 0.0206826
beta3_td = 0.0211664
beta4_td = 0.0200515
beta5_td = -0.0230394
alpha_td = -0.434975

...
...
...

~~~
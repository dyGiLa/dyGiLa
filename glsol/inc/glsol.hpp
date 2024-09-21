#ifndef GLSOL_HPP
#define GLSOL_HPP

#define USE_PARIO 
#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
//#include <math.h>
#include <assert.h>

#include "plumbing/hila.h"
#include "plumbing/fft.h"

#include "matep.hpp"
#include "dyGiLa_config.hpp"

// Definition of the field that we will use
using real_t = float;                          // or double ?
using phi_t = Matrix<3,3,Complex<real_t>>;     // saves the trouble of writing this every time

// Container for simulation parameters and methods
class glsol{

public:
  glsol() = default;                     // default constructor
  
  // read configration file and initiate scaling_sim.config 
  const std::string allocate(const std::string &fname, int argc, char **argv);
  
  void initialize();
  
  void initializeT();
  void initializep();
  
  void point_params(real_t T, real_t p, real_t beta[6]);
  
  void write_moduli();
  void write_energies();
  void write_positions();
  void write_phases();
  
  void next();
  void next_bath();
  void next_bath_UniT_quench();  
  void nextT();
  void hotbloob();
  
  Field<phi_t> A;
  Field<phi_t> pi;

  Field<real_t> T;
  Field<real_t> dT;
  Field<real_t> dT_from_local_TAB;  
  Field<real_t> p;
  
  real_t t;
  real_t tc = 0;

  Matep MP;

  dyGiLaConf config;

  std::vector<real_t> t_v;
  std::vector<real_t> T_v;
  std::vector<real_t> p_v;
    
};

#endif

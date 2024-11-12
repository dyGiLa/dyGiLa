#ifndef GLSOL_HPP
#define GLSOL_HPP

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
  void initializeH();
  
  void point_params(real_t T, real_t p, real_t beta[6]);
  
  void write_moduli();
  void write_energies();
  void write_positions();
  void write_phases();
  
  //void next();
  void next_bath();
  //void next_bath_UniT_quench();  
  void nextT();
  void hotbloob();

  phi_t shift(phi_t phip ,phi_t phi0, phi_t phim, CoordinateVector p, Direction d, int dir);  
  phi_t shift_talk(phi_t phip ,phi_t phi0, phi_t phim, CoordinateVector p, Direction d, int dir);
  phi_t periodic(phi_t phi1);
  phi_t maximal();
  phi_t BphaseBoundary();
  phi_t AphaseBoundary();
  phi_t normalBoundary();
  phi_t RobinBoundary(phi_t phiI, int oorN, real_t bt, int face);

  real_t shiftT(real_t Tp, real_t T0, real_t Tm, CoordinateVector p, Direction d, int dir);
  real_t periodic_T(real_t T1);
  real_t fixedT();
  
  Field<phi_t> A;
  Field<phi_t> pi;

#ifdef T_FIELD
  Field<real_t> T;
  Field<real_t> dT;
#else
  real_t T;
#endif

  real_t p;
  // Field<real_t> dT_from_local_TAB;  //check if it is useful or not
  //Field<real_t> p; // it can be removed now but it can be a work for the future

  real_t H[3]; //magnetic field;
  real_t gh=0.852;
  real_t gz=0.022;

  
  real_t t;
  real_t tc = 0;

  Matep MP;

  dyGiLaConf config;

  std::vector<real_t> t_v;
  std::vector<real_t> T_v;
  std::vector<real_t> p_v;
    
};

#endif

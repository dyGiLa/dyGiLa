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


// Define convenient enum for addressing the reductions - could use numbers too

/* enumeration counter type  for matrix reduction */
// enum class Matreduc { i_sumA, N_MatREDUCTION };
enum class matreduc { i_sumA, N_MatREDUCTION };

/* enumeration counter type for scalar reduction */
//enum class reduc {
enum {
  i_sumgapA,
  i_suma, i_sumb1, i_sumb2, i_sumb3, i_sumb4, i_sumb5,
  i_suma_we, i_sumb1_we, i_sumb2_we, i_sumb3_we, i_sumb4_we, i_sumb5_we,
  i_sumk1, i_sumk2, i_sumk3, i_sumk1_we, i_sumk2_we, i_sumk3_we,
  i_sumkin, i_sumkin_we,
  N_REDUCTION
};  

// Container for simulation parameters and methods
class glsol{

public:
  glsol() = default;                     // default constructor
  
  // read configration file and initiate scaling_sim.config 
  const std::string allocate(const std::string &fname, int argc, char **argv);

  // OP field initialization
  void initialize();

  // T-field initialization  
  void initializeT();

  // p-field initialization
  // void initializep();
  
  // H-field initialization    
  void initializeH();

  // bulk coefficients calculating tool
  void point_params(real_t T, real_t p, real_t beta[6]);

  void write_moduli();
  void write_energies();
  void write_positions();
  void write_phases();
  
  void phaseMarking();
  
  void next();
  void next_bath();
  void next_bath_UniT_quench();
  void next_bath_UniT_quench_Hfield();
  void next_bath_hotblob_quench_Hfield();      
  //void nextT();
  
  Field<phi_t> A;
  Field<phi_t> pi;

  Field<real_t> T;
  Field<real_t> dT;
  // Field<real_t> dT_from_local_TAB;  
  // Field<real_t> p;

  Field<Vector<3,real_t>> H; // H-field, 3-component column vector field

  Field<real_t> phaseMarker;
  
  real_t t;
  real_t tc = 0;

  matep::Matep MP;

  dyGiLaConf config;

  std::vector<real_t> t_v;
  std::vector<real_t> T_v;
  std::vector<real_t> p_v;
    
};

#endif

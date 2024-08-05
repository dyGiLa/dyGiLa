#ifndef GLSOL_HPP
#define GLSOL_HPP

#define _USE_MATH_DEFINES
#define USE_ASCENT 
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
  void nextT();
  void hotbloob();
  
  Field<phi_t> A;
  Field<phi_t> pi;

  Field<real_t> T;
  Field<real_t> dT;
  Field<real_t> p;
  
  real_t t;
  real_t tc = 0;

  std::vector<real_t> t_v;
  std::vector<real_t> T_v;
  std::vector<real_t> p_v;
  
  struct config {
      int lx;
      int ly;
      int lz;
      real_t dx;
      real_t dt;
      real_t dtdxRatio;
      
      real_t tStart;
      real_t tEnd;

      real_t tdif;
      real_t difFac;
      real_t tdis;

      //real_t gamma;
      Complex<real_t> gamma;
      Complex<real_t> gamma1;
      Complex<real_t> gamma2;
      int gammaoffc;
      
      int initialCondition;
      real_t variance_sigma;
      
      int seed;
      real_t IniMod;
      real_t Inilc;

      int item;
    

    /* These are new, need look in... */
      int initialConditionT;
      real_t IniT;
      real_t ampT;
      real_t sigTx;
      real_t sigTy;
      real_t sigTz;
      int initialConditionp;
      real_t Inip;
    /* These are new, need look in... */
    
    
      real_t T;
      real_t dT_from_TAB;
      real_t p;
    
      real_t alpha;
      real_t beta1;
      real_t beta2;
      real_t beta3;
      real_t beta4;
      real_t beta5;

      real_t tStats;
      real_t nOutputs;

      std::fstream stream;
      
      // xdmf file fstream
      std::fstream xdmf_out;
      std::fstream xml_out;      
      std::string  xmf2_fname;

      /*----------------------------------------*/
      /*          switches declearition         */
      /*----------------------------------------*/
      //unsigned int A_matrix_output;
      /*----------------------------------------*/
      /*        switches declearation end       */
      /*----------------------------------------*/      

      int positions;
      int npositionout;
      
      int boundaryConditions;
      int BCs1;
      int BCs2;
      // int Wn;
      int BCchangec;
      
      // real_t BLeft_11,BLeft_22,BLeft_33,
      // 	     BRight_11,BRight_22,BRight_33;

      int useTbath;

      /*----------------------------------------*/
      /*       insitu parameter declearition    */
      /*----------------------------------------*/
      // unsigned int do_gapA_clip,
      //              do_gapA_isosurface,
      //              do_gapA_3slice,
      //              do_fe_slice,
      //              do_gapA_slice;

      // unsigned int hdf5_A_matrix_output,
      //              hdf5_trA_output,
      //              hdf5_eigvA_output,
      //              hdf5_mass_current_output,
      //              hdf5_spin_current_output;
	
      
      real_t clamp_bias_gapMin, clamp_bias_gapMax;
      // real_t clamp_fed_Min, clamp_fed_Max;      
      real_t camera_azi, camera_ele;
      /*----------------------------------------*/
      /*       insitu parameter end             */
      /*----------------------------------------*/
      
      int write_phases;
      int write_eigen;

    /* These are new, need look in... */
      int evolveT;
      int Tevolvetype;
      real_t startdiffT;
      real_t diffT;
      int bloob_after;
      real_t theat;
    /* These are new, need look in... */    
    
  } config;
  
};

#endif

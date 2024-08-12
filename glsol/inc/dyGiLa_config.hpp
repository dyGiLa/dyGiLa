#ifndef DYGILA_CONFIG_HPP
#define DYGILA_CONFIG_HPP

#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

using real_t = float;                     


// parameters configuraion class of dyGiLa
struct dyGiLaConf {
  
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
        
      // real_t T;
      // real_t dT_from_TAB;
      // real_t p;
    
      real_t alpha;
      real_t beta1;
      real_t beta2;
      real_t beta3;
      real_t beta4;
      real_t beta5;

      int initialConditionT;
      real_t IniT;
      real_t ampT;
      real_t sigTx;
      real_t sigTy;
      real_t sigTz;
      int initialConditionp;
      real_t Inip;
 
      real_t tStats;
      real_t nOutputs;

      std::fstream stream;
      
      // xdmf file fstream
      std::fstream xdmf_out;
      std::fstream xml_out;      
      std::string  xmf2_fname;

      int positions;
      int npositionout;
      
      int boundaryConditions;
      int BCs1;
      int BCs2;
      // int Wn;
      int BCchangec;
      
      int useTbath;

      int write_phases;
      int write_eigen;

      int evolveT;
      int Tevolvetype;
      real_t startdiffT;
      real_t diffT;
      int bloob_after;
      real_t theat;
    
      /*----------------------------------------*/
      /*     parallel IO control parameters     */
      /*----------------------------------------*/
      unsigned int A_matrix_output;
      // unsigned int do_gapA_clip,
      //              do_gapA_isosurface,
      //              do_gapA_3slice,
      //              do_fe_slice,
      //              do_gapA_slice;

      unsigned int hdf5_A_matrix_output,
      //              hdf5_trA_output,
      //              hdf5_eigvA_output,
                   hdf5_mass_current_output,
                   hdf5_spin_current_output;
	
      
      real_t clamp_bias_gapMin, clamp_bias_gapMax;
      real_t clamp_fed_Min, clamp_fed_Max;      
      real_t camera_azi, camera_ele;
      /*----------------------------------------*/
      /*    parallel IO parameter end           */
      /*----------------------------------------*/
};
  

#endif

#ifndef DYGILA_CONFIG_HPP
#define DYGILA_CONFIG_HPP

#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "plumbing/hila.h"

using real_t = float;                     

// parameters configuraion class of dyGiLa
struct dyGiLaConf {
  
      int lx;
      int ly;
      int lz;
      real_t dx;
      real_t dt;
      real_t dtdxRatio;

      // tThermalizationWaiting, initilial configuration thermalizing time
      // quench time, measured from Tc to 0
      // Ttd_Qend, detinated temperature Ttd_Qend = T_Qend/Tc
      // has1stQStop, has 1st quench stop point
      // Ttd_Q1st, 1st quench stop Temp, Ttd_Q1st = T_Q1st/Tc
      // tQ1Waiting, waiting time of 1st quench stop Temp
      real_t tThermalizationWaiting;
      real_t tauQ1;
      real_t tauQ2;  
      real_t tauQ;
      unsigned int has1stQStop;
      real_t Ttd_Q1st;
      real_t tQ1Waiting;
      real_t Ttd_Qend;
  
      real_t tStart;
      real_t tEnd;

      real_t tdif;
      real_t difFac;
      real_t tdis;

      //real_t gamma;
      unsigned int TDependnetgamma;
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

      // T-filed initialing & control parameters
      int initialConditionT;
      real_t IniT;
      real_t ampT; //sine wave profile
      /*-----------------------------------------*/
      /* hot blob Temeprature profile parameters */
      /* Ttdbx mean T in unit of Tc              */
      real_t Ttdb1;
      real_t Ttdb0;
      real_t t1;
      /*-----------------------------------------*/  
 

      // p-filed initialing & control parameters
      int initialConditionp;
      real_t Inip;

      // H-filed initialing & control parameters
      unsigned int withHfield;  
      int initialConditionH;
      Vector<3, real_t> InitH;
  
      real_t tStats;
      real_t nOutputs;

      std::fstream stream;
      
      std::string  xmf2_fname;

      int positions;
      int npositionout;
      
      int boundaryConditions;
      int BCs1;
      int BCs2;
      // int Wn;
      int BCchangec;
      
      int useTbath;
      real_t Tbath_start;

      int write_phases;
      int write_eigen;

      real_t ptol;
  
      int evolveT;
      int Tevolvetype;
      real_t startdiffT;
      real_t diffT;
      int bloob_after;
      real_t theat;
    
      /*----------------------------------------*/
      /*     parallel IO control parameters     */
      /*----------------------------------------*/
      unsigned int do_gapA_clip;
      real_t gapA_clip1_point_x, gapA_clip1_point_y, gapA_clip1_point_z,
             gapA_clip1_norm_x, gapA_clip1_norm_y, gapA_clip1_norm_z;
      real_t gapA_clip2_point_x, gapA_clip2_point_y, gapA_clip2_point_z,
             gapA_clip2_norm_x, gapA_clip2_norm_y, gapA_clip2_norm_z;

      unsigned int do_gapA_slice;
      real_t gapA_slice1_point_x, gapA_slice1_point_y, gapA_slice1_point_z,
             gapA_slice1_norm_x, gapA_slice1_norm_y, gapA_slice1_norm_z;
      real_t gapA_slice2_point_x, gapA_slice2_point_y, gapA_slice2_point_z,
             gapA_slice2_norm_x, gapA_slice2_norm_y, gapA_slice2_norm_z;
 
  
      unsigned int do_fed_clip;
      real_t fed_clip_point_x, fed_clip_point_y, fed_clip_point_z,
             fed_clip_norm_x, fed_clip_norm_y, fed_clip_norm_z;

      unsigned int do_Temperature_clip;
      real_t Temperature_clip_point_x, Temperature_clip_point_y, Temperature_clip_point_z,
	     Temperature_clip_norm_x, Temperature_clip_norm_y, Temperature_clip_norm_z,
             Temperature_clamp;

      unsigned int do_Temperature_slice;
      real_t Temperature_slice_point_x, Temperature_slice_point_y, Temperature_slice_point_z,
	     Temperature_slice_norm_x, Temperature_slice_norm_y, Temperature_slice_norm_z;
  
  
      unsigned int do_Temperature_isosurface;
      std::vector<double> Temperature_iso_values_vector;

      unsigned int do_gapA_isosurface;
      std::vector<double> iso_values_vector;

      unsigned int do_phaseMarker_slice;
      real_t pMarker_slice_point_x, pMarker_slice_point_y, pMarker_slice_point_z,
	     pMarker_slice_norm_x, pMarker_slice_norm_y, pMarker_slice_norm_z;
    
      //              do_gapA_3slice,
      //              do_fe_slice,
      //              do_gapA_slice;

      unsigned int hdf5_A_matrix_output;
      real_t hdf5Ststart, hdf5Stend;	
      //              hdf5_trA_output,
      //              hdf5_eigvA_output,
     unsigned int  hdf5_mass_current_output,
                   hdf5_spin_current_output;
	
      
      real_t clamp_bias_gapMin, clamp_bias_gapMax;
      real_t clamp_bias_fed_Min, clamp_bias_fed_Max;      
      real_t camera1_azi, camera1_ele,
             camera2_azi, camera2_ele;
      /*----------------------------------------*/
      /*    parallel IO parameter end           */
      /*----------------------------------------*/
};
  

#endif

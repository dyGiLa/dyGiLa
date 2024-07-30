#ifndef PARIO_HPP
#define PARIO_HPP

#define _USE_MATH_DEFINES
//#define USE_ASCENT 
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
#include "glsol.hpp"

/*-----------------------------------------------------------------------*/
/***      include Ascent & Canduit for in situ rank rendering        *****/
/*-----------------------------------------------------------------------*/
#if defined USE_ASCENT

#include "ascent.hpp"
#include "conduit_blueprint.hpp"

#endif
/*-----------------------------------------------------------------------*/
/***               include Ascent & Canduit end here                 *****/
/*-----------------------------------------------------------------------*/

// Definition of the field that we will use
using real_t = float;                          // or double ?
// using phi_t = Matrix<3,3,Complex<real_t>>;     // saves the trouble of writing this every time

// class tamplate for handling parallel stream
// template paramter is the dimension of order paramter matrix e.g., 1, 2, 3
template <unsigned int OPDim> 
class parIO{

public:
  parIO() = default;                     // default constructor
  
  Matep matep;
  
  /*----------------------------------------*/
  /*  parallel-IO in-situ via Ascent        */
  /*----------------------------------------*/

  // called in main() before t-while-loop started
  void insitu_hdf5xdmf(glsol &);

  // called in main()  
  void insitu_initialize(glsol &);
  void insitu_execute(glsol &);
  void insitu_close();  

  // called in insitu_initialize()
  void insitu_createMesh(glsol &);
  void insitu_defineActions(glsol &);

  /*----- fields declearations -----*/
  
  Field<real_t> gapA;
    // Field<real_t> feDensity;
    // Field<real_t> trA_re, trA_im;
    // Field<real_t> u11, u12, u13, u21, u22, u23, u31, u32, u33;
    // Field<real_t> v11, v12, v13, v21, v22, v23, v31, v32, v33;
    // Field<real_t> eigAv1, eigAv2, eigAv3;
    // Field<real_t> jm1, jm2, jm3;
    // Field<real_t> phaseExpModulus, phaseExpAngle,/*acosphi, asinphi,*/ phaseExp2Re, phaseExp2Im;
    // Field<real_t> js11, js21, js31,
    //               js12, js22, js32,
    //               js13, js23, js33;
   
  std::vector<real_t> gapAOrdered;
    // std::vector<real_t> feDensityOrdered;
    // std::vector<real_t> trA_reOrdered, trA_imOrdered;
    // std::vector<real_t> u11Ordered, u12Ordered, u13Ordered,
    //                     u21Ordered, u22Ordered, u23Ordered,
    //                     u31Ordered, u32Ordered, u33Ordered;
    // std::vector<real_t> v11Ordered, v12Ordered, v13Ordered,
    //                     v21Ordered, v22Ordered, v23Ordered,
    //                     v31Ordered, v32Ordered, v33Ordered;
    // std::vector<real_t> eigAv1Ordered, eigAv2Ordered, eigAv3Ordered;
    // std::vector<real_t> jm1Ordered, jm2Ordered, jm3Ordered;
    // std::vector<real_t> phaseExpModulusO, phaseExpAngleO,/*acosphiO, asinphiO,*/ phaseExp2ReO, phaseExp2ImO;
    // std::vector<real_t> js11O, js21O, js31O,
    //                     js12O, js22O, js32O,
    //                     js13O, js23O, js33O;
  
    /*--------------------------------*/
  
  long long ghostVolume;
  long long latticeVolumeWithGhost;
  long long latticeVolume;
  
  unsigned char *ghostCellsMask;

  ascent::Ascent insitu;
  conduit::Node ascent_options;
  conduit::Node actions;
  conduit::Node mesh;

  /*----------------------------------------*/
  /*     In-situ declearations end here     */
  /*----------------------------------------*/
  
};

#endif

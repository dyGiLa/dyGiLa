#ifndef PARIO_HPP
#define PARIO_HPP

#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>

#include "plumbing/hila.h"
#include "plumbing/fft.h"

#include "matep.hpp"
#include "glsol.hpp"

#include "ascent.hpp"
#include "conduit_blueprint.hpp"


// Definition of the field that we will use
using real_t = float;                          // or double ?
// using phi_t = Matrix<3,3,Complex<real_t>>;     // saves the trouble of writing this every time

class parIO{

public:
  parIO() = default;                     // default constructor
  
  matep::Matep matep;
  
  /*----------------------------------------*/
  /*  parallel-IO memember functions        */ 
  /*----------------------------------------*/

  // called in main() before t-while-loop started
  void xdmf(glsol &);
  void xml(glsol &);

  // called in main()  
  void init(glsol &);
  void pstream(glsol &);
  void shutdown();

  // xdmf file fstream
  std::fstream xdmf_out;
  std::fstream xml_out;

private:
  
  /* functions called pario engine_init() */
  /*  mesh descriptions and definations   */
  void describeMesh(glsol &);
  void describeMesh_addGhost_verify();
  
  void describeMesh_gapA_FEDensity();
  void describeMesh_Temperature();
  void describeMesh_phaseMarker();
  void describeMesh_massCurrent();
  void describeMesh_spinCurrent();  
  void describeMesh_AMatrix();  

  /* actions definations */
  void defineActions_insitu(glsol &);

  void defineActions_massCurrent();
  void defineActions_spinCurrent();
  void defineActions_AMatrix();    
  
  void defineActions_printTree();  

  /*----- fields declearations -----*/
  
  Field<real_t> gapA;
  Field<real_t> feDensity;
    // Field<real_t> trA_re, trA_im;
  Field<real_t> u11, u12, u13, u21, u22, u23, u31, u32, u33;
  Field<real_t> v11, v12, v13, v21, v22, v23, v31, v32, v33;
    // Field<real_t> eigAv1, eigAv2, eigAv3;
  Field<real_t> jm1, jm2, jm3;
  Field<real_t> phaseExpModulus, phaseExpAngle/*acosphi, asinphi,*/
                ,phaseExp2Re, phaseExp2Im;
  Field<real_t> js11, js21, js31,
                js12, js22, js32,
                js13, js23, js33;
   
  std::vector<real_t> gapAContainer;
  std::vector<real_t> feDensityContainer;
  std::vector<real_t> Temperature;
  std::vector<real_t> phaseMarker;  

  // std::vector<real_t> trA_reContainer, trA_imContainer;
  std::vector<real_t> u11Container, u12Container, u13Container,
                      u21Container, u22Container, u23Container,
                      u31Container, u32Container, u33Container;
  std::vector<real_t> v11Container, v12Container, v13Container,
                      v21Container, v22Container, v23Container,
                      v31Container, v32Container, v33Container;
    // std::vector<real_t> eigAv1Container, eigAv2Container, eigAv3Container;
  std::vector<real_t> jm1Container, jm2Container, jm3Container;
  std::vector<real_t> phaseExpModulusContainer, phaseExpAngleContainer /*acosphiContainer, asinphiContainer,*/
                      ,phaseExp2ReContainer, phaseExp2ImContainer;
  std::vector<real_t> js11Container, js21Container, js31Container,
                      js12Container, js22Container, js32Container,
                      js13Container, js23Container, js33Container;
  
  /*--------------------------------*/
  
  long long ghostVolume;
  long long latticeVolumeWithGhost;
  long long latticeVolume;
  
  unsigned char *ghostCellsMask;

  ascent::Ascent pio;
  conduit::Node pio_options;
  conduit::Node actions;
  conduit::Node mesh;
  
  /*----------------------------------------*/
  /*      pario declearations end here      */
  /*----------------------------------------*/
  
};

#endif

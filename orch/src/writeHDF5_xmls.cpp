#define USE_PARIO 
#define USE_MPI 
// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <string>
// #include <assert.h>

//#include "plumbing/hila.h"
//#include "plumbing/globals.h" 

#include "glsol.hpp"
//#include "matep_namespace_utils.hpp"
#include "orch.hpp"

#if defined USE_PARIO 
#include "pario.hpp"
#endif

namespace orch {
  
void writeHDF5_xmls(glsol &gl, parIO &paraio) {
  //xml files for MetaData.
  if (
      (gl.config.hdf5_A_matrix_output        == 1)
      || (gl.config.hdf5_mass_current_output == 1)
      || (gl.config.hdf5_spin_current_output == 1)
      || (gl.config.hdf5_pMarker_output == 1)
      ) {
    if (gl.config.hdf5_A_matrix_output == 1) paraio.xml_Amatrix(gl);
    if (gl.config.hdf5_pMarker_output == 1) paraio.xml_pMarker(gl);
    if (gl.config.hdf5_mass_current_output == 1) paraio.xml_massCurrent(gl);
    if (gl.config.hdf5_spin_current_output == 1) paraio.xml_spinCurrent(gl); 
  }

} // writeHDF5_xmls() func endds here

} // orch namespace ends here

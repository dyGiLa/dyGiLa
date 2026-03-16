//#define USE_PMD_GAMMA
//#define USE_PARIO 
#define USE_MPI 
// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <string>
// #include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/globals.h" 

#include "glsol.hpp"
//#include "matep_namespace_utils.hpp"
//#include "dyGiLa.hpp"
#include "dyGiLa_utils.hpp"

// #if defined USE_PARIO 
// #include "pario.hpp"
// #endif

namespace dyGiLa {
  namespace utils {
    const bool latticeSaturationCheck(glsol &gl)
       {
	 const bool saturated = ((gl.pM5VolR >= 0.9999) || (gl.pM9VolR >= 0.9999)) ? true : false;
	 if (saturated == true)
	   {
	     hila::out0 << "\n-----------------------------------------------------------------------------"
	                << "\n"
	                << "lattice is saturated, terminate iteration! "
	                << "\n-----------------------------------------------------------------------------"
			<< std::endl;
	   }
	 return saturated;
       } // latticeSaturationCheck() func ends here
  } // dyGiLa::utils:: namespace ends here
} // dyGiLa namespace ends here

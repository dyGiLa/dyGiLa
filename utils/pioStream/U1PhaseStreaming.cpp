#define USE_BOTHSIDE_GHOSTS
#define USE_ADGRZ
// #define USE_MPI 
// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <string>
// #include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/fft.h"
//#include "plumbing/memalloc.h"

#include "glsol.hpp"
//#include "matep.hpp"
#include "pario.hpp"

// #include "ascent.hpp"
// #include "conduit_blueprint.hpp"

void parIO::U1PhaseStreaming(glsol &sol) {

  onsites(ALL)
    {
      Complex<real_t> detAwT_X = sol.AwT[X].det_laplace();

      /* -----------------------------------------------------------
       * det(AwT) = e^{3\phi} gap^3/3sqrt(3) for standard B-phase OP,
       * for phases breaking TR-symmetry, such as A-phase, det(AwT) = 0 (or << O(1)),
       * before introducing A-phase U(1) computation, U1_3phi is simply given as zero.
       * -----------------------------------------------------------
       */
      U1_3phi[X] = (detAwT_X.abs() >= 1e-3)
	           ?(detAwT_X/detAwT_X.abs()).arg()
	           : 0.;
    }

  U1_3phi.copy_local_data_with_halo(U1_3phiContainer);
} // U1PhaseStreaming(glsol &) end here


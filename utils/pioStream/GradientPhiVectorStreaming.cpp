//#define USE_BOTHSIDE_GHOSTS
//#define USE_ADGRZ
//#define PI 3.141592653589793
// #define USE_MPI 
// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <string>
// #include <assert.h>

#include "plumbing/hila.h"
#include "plumbing/coordinates.h"
//#include "plumbing/fft.h"
//#include "plumbing/memalloc.h"

#include "glsol.hpp"
#include "matep.hpp"
#include "pario.hpp"

// #include "ascent.hpp"
// #include "conduit_blueprint.hpp"

void parIO::GradientPhiVectorStreaming(glsol &sol) {
  onsites(ALL)
    {
      matep::Matep MPonsites;

      Vector<3,Complex<real_t>> ApjAd;
      const real_t invDelta2 = 1.f/(((sol.AwT[X]*(sol.AwT[X].dagger())).trace()).real());

      foralldir(j) {
	ApjAd[j] = (1.f/(2.*sol.config.dx)) * ((sol.AwT[X] * sol.AwT[X + j].dagger()).trace()
			                       - (sol.AwT[X] * sol.AwT[X - j].dagger()).trace());
	
      } // ApjAd for j = 0,1,2

      GPhi_1[X] = -invDelta2 * ApjAd[0].imag();
      GPhi_2[X] = -invDelta2 * ApjAd[1].imag();
      GPhi_3[X] = -invDelta2 * ApjAd[2].imag();      
                          
    } // onsite() block done

} // lVectorStreaming(glsol &) end here


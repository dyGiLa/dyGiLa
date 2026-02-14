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
//#include "plumbing/coordinates.h"
//#include "plumbing/fft.h"
//#include "plumbing/memalloc.h"

#include "glsol.hpp"
#include "matep.hpp"
#include "pario.hpp"

// #include "ascent.hpp"
// #include "conduit_blueprint.hpp"

void parIO::lVectorStreaming(glsol &sol) {
  onsites(ALL)
    {
      matep::Matep MPonsites;

      const Matrix<3,3,real_t> ImAdA_X = (sol.AwT[X].dagger() * sol.AwT[X]).imag();
      const real_t invDelta2 = 1.f/(((sol.AwT[X]*(sol.AwT[X].dagger())).trace()).real());

      real_t l1, l2, l3;

      l1 = invDelta2 * MPonsites.epsilon(1,2,0) * ImAdA_X.e(1,2)
	   + invDelta2 * MPonsites.epsilon(2,1,0) * ImAdA_X.e(2,1);

      l2 = invDelta2 * MPonsites.epsilon(0,2,1) * ImAdA_X.e(0,2)
	   + invDelta2 * MPonsites.epsilon(2,0,1) * ImAdA_X.e(2,0);

      l3 = invDelta2 * MPonsites.epsilon(1,0,2) * ImAdA_X.e(1,0)
	   + invDelta2 * MPonsites.epsilon(0,1,2) * ImAdA_X.e(0,1);
            
      /* l_1 = \epsilon_ij0 Im(Ad A)_ij */
      //foralldir(j) foralldir(i) { l1 += /*invDelta2 * MPonsites.epsilon(i,j,0) */ ImAdA_X.e(i,j); }

      /* l_2 = \epsilon_ij1 Im(Ad A)_ij */
      // foralldir(j) foralldir(i) { l2 += /*invDelta2 */ MPonsites.epsilon(i,j,1) * ImAdA_X.e(i,j); }

      // /* l_3 = \epsilon_ij2 Im(Ad A)_ij */
      // foralldir(j) foralldir(i) { l3 += invDelta2 * MPonsites.epsilon(i,j,2) * ImAdA_X.e(i,j); }

      const real_t l_Squ = l1*l1 + l2*l2 + l3*l3;
      // const real_t l_Squ = invDelta2 * (MPonsites.epsilon(0,1,2) * ImAdA_X.e(0,1) + MPonsites.epsilon(1,0,2) * ImAdA_X.e(1,0));      

      if (
	  ( l_Squ <= 1.f + sol.config.lVec_SqlTol )
	  //&& ( l_Squ > 1.f - sol.config.lVec_SqlTol )
	 )
	{
	  lsq[X] = l_Squ; //l_1[X] = l1; l_2[X] = l2; l_3[X] = l3;
          if (sol.config.hdf5_lVector_output == 1)
	    {l_1[X] = l1; l_2[X] = l2; l_3[X] = l3;}	  
	}
      else
	{
	  lsq[X] = 0.f; //l_1[X] = 0.f; l_2[X] = 0.f; l_3[X] = 0.f;
	  if (sol.config.hdf5_lVector_output == 1)
	    {l_1[X] = 0.f; l_2[X] = 0.f; l_3[X] = 0.f;}	  
	}
                     
    } // onsite() block done

} // lVectorStreaming(glsol &) end here


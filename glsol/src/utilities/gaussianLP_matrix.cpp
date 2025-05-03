#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>

#include "plumbing/hila.h"
#include "plumbing/fft.h"

#include "glsol.hpp"
#include "matep.hpp"


void glsol::GaussianLPfilter_matrix(Field<phi_t> &AwT /* Weierstrass Transformed */) {

     /* --------------------------------------------------------------------------------------------
      * declared and defined Gaussian Low Pass filter member in glsol.hpp.
      * This specific implementation is for Filed<phi_t> type, 
      * which is 3x3 complex valued matrix field.
      * 
      * The algorithm is basded on fact that Gaussian Low Pass filter 
      * is Weiestraas transformation (WT), which suppress 
      * shortwave length noise expotentially. 
      * Weiestrass transformation of a given spatial configuration (field containing shortwave noise) 
      * is solution of an equvlent heat equation. Therefore, algorithm used in this filter 
      * constructs rough solution of this equvlent heat equtation by 
      * doing few step of finite difference iteration.
      *
      * Iterations in this func are adding correction to last step Field<phi_t> object.
      * correction is constrcted with finite difference in symmetry fishion i.e.,
      * Ax^(m+1) = Ax^(m) + c2 (Ax-1^(m) -2Ax^(m) + Ax+1^(m)) along every direction.
      * In 3D case, the coefficient is 3 x 2c2 = 6c2.
      *
      * Refine above formulism, we have the algorithm:
      * Ax^(m+1) = (Ax^(m) - 2*c2 Ax^(m)) + c2*(Ax-1^(m) + Ax+1^(m)) along every direction.
      *
      * Then one has c1 = 1-6c2, and c2 is the inversed squared time step-length.
      * such asc1 = 1/3, c2 = 1/9.
      * ---------------------------------------------------------------------------------------------
      */

     for (unsigned int iter=0; iter < config.numIterGLPfilter; iter++)
       {
        for (Parity par : {EVEN,ODD})
	  {
	   onsites(par)
	     {
              AwT[X] *= config.GLPfc1;
	      foralldir(d) AwT[X] += config.GLPfc2 * (AwT[X-d] + AwT[X+d]);
	     }
	  }
       
       // onsites(ALL) {
       //    As2[X] = config.smearC2 * As1[X];
       //    foralldir(d) As2[X] += config.smearC2 * (As1[X-d] + As1[X+d]);
       // }
       // hila::swap(As1,As2);   // flips content without copying it
       } // AwT is filtered, shortwave noise suppressed  - measure phase from it

  
} // GaussianLPfilter_matrix() function ends here


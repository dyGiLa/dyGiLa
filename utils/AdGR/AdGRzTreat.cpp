#define USE_MPI 
//#include <sstream>
//#include <iostream>
//#include <iomanip>
//#include <fstream>
//#include <string>
//#include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/fft.h"
#include "plumbing/coordinates.h"

#include "glsol.hpp"
//#include "matep.hpp"


void glsol::AdGRzTreat(phi_t &A_Xm, phi_t &A_Xp, const real_t &abt_ratio) {

    if ( X.coordinate(e_z) == 0 )
      {
        const real_t trcoef=(1.-abt_ratio)/(1.+abt_ratio);
	foralldir(col)
	  {
           if(col == 2)
	     { foralldir(row){ A_Xm.e(row, col)=0.0; } }
	   else
	     { foralldir(row){ A_Xm.e(row, col)=A[X+e_z].e(row, col)*trcoef;} }
	  } // col loop ends here
      }
    else if ( X.coordinate(e_z) == (config.lz-1) )
      {
        const real_t trcoef=(1.+abt_ratio)/(1.-abt_ratio);
	foralldir(col)
	  {
           if(col == 2)
	     { foralldir(row){ A_Xp.e(row, col)=0.0; } }
	   else
	     { foralldir(row){ A_Xp.e(row, col)=A[X-e_z].e(row, col)*trcoef;} }
	  } // col loop ends here
      }
   
} // AdGRzTreat() ends here


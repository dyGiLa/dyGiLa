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

#include "glsol.hpp"
#include "matep.hpp"

void glsol::nextT() {

  Field<real_t> deltaT;
  real_t alpha=0.1;
  real_t heatfactor=3.0;

  real_t time_steps=1.0;
  real_t dt=config.dt;

#ifdef T_FIELD
  
  switch (config.Tevolvetype){
   case 0: {
     // case heat
     dt=config.dt/heatfactor;
     time_steps=heatfactor;
    
     break;
   }
   case 1: {

     // case wave
     dt=config.dt;
     time_steps=1;
    
     break;
   }
  } // config.Tevolvetype block
  
  for(int i=0; i<time_steps; i++)
    {
  
      onsites (ALL) {

	T[X] += dt * dT[X];
      }

      onsites (ALL) {

	//real_t diff2 = alpha * (1.0/(config.dx*config.dx)) * (T[X + e_x] + T[X - e_x]
	//						     + T[X + e_y] + T[X - e_y]
	//						     + T[X + e_z] + T[X - e_z]
	//						     - 6.0*T[X]); 

	real_t diff2 = alpha * (1.0/(config.dx*config.dx)) * (shiftT(T[X+e_x],T[X],T[X-e_x],X.coordinates(),e_x,1) + shiftT(T[X+e_x],T[X],T[X-e_x],X.coordinates(),e_x,0)
							   + shiftT(T[X+e_y],T[X],T[X-e_y],X.coordinates(),e_y,1) + shiftT(T[X+e_y],T[X],T[X-e_y],X.coordinates(),e_y,0)
							   + shiftT(T[X+e_z],T[X],T[X-e_z],X.coordinates(),e_z,1) + shiftT(T[X+e_z],T[X],T[X-e_z],X.coordinates(),e_z,0)
							   - 6.0*T[X]);

	
	switch (config.Tevolvetype){
	case 0: {
	  
	  dT[X] = diff2;
	
	  break;
	}
	case 1: {
	  
	  dT[X] = dT[X] + diff2 * dt;

	  break;
	}
	}
      } // onsite(ALL) block

    } // for loop block
#else



  
#endif  
} // nextT() ends here


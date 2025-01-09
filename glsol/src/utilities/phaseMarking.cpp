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


void glsol::phaseMarking() {

  onsites (ALL) {

    real_t R1{0},R2{0},R3{0},R4{0},R5{0};
    
    // reduced OP matrix; TrAr.Ar^dagger = 1
    phi_t Ar=A[X]/sqrt((A[X]*A[X].dagger()).trace());

    R1 = ((Ar*Ar.transpose()).trace()).squarenorm();

    R2 = real(((Ar*Ar.dagger()).trace()*(Ar*Ar.dagger()).trace()));

    R3 = real(((Ar*Ar.transpose()*Ar.conj()*Ar.dagger()).trace()));

    R4 = real(((Ar*Ar.dagger()*Ar*Ar.dagger()).trace()));

    R5 = real(((Ar*Ar.dagger()*Ar.conj()*Ar.transpose()).trace()));

    /* >>>>>>>>> phase marking logic <<<<<<<<  */
    if (
	(abs(R1-0.0) <= config.ptol)
	&& (abs(R3-0.0) <= config.ptol)
       )
      {
	// A-phase OP onsite
	if (
	    (abs(R2-1.0) <= config.ptol)
	    && (abs(R4-1.0) <= config.ptol)
	    && (abs(R5-1.0) <= config.ptol)
	   )
	  phaseMarker[X] = 4.0f;
	else
	  phaseMarker[X] = 0.0f;	  
      }
    else if (abs(R1-1.0) <= config.ptol)
      {
	// B-phase OP onsite
        if (
            (abs(R2-1.0) <= config.ptol)
	    && (abs(R3-(1.0/3.0)) <= config.ptol)
	    && (abs(R4-(1.0/3.0)) <= config.ptol)
	    && (abs(R5-(1.0/3.0)) <= config.ptol)
           )
	  phaseMarker[X] = 2.0f;
	// planar-phase OP onsite
	else if (
                 (abs(R2-1.0) <= config.ptol)
	         && (abs(R3-(1.0/2.0)) <= config.ptol)
	         && (abs(R4-(1.0/2.0)) <= config.ptol)
	         && (abs(R5-(1.0/2.0)) <= config.ptol)
                )
	  phaseMarker[X] = 1.0f;
	// polar-phase OP onsite
	else if (
                 (abs(R2-1.0) <= config.ptol)
	         && (abs(R3-1.0) <= config.ptol)
	         && (abs(R4-1.0) <= config.ptol)
	         && (abs(R5-1.0) <= config.ptol)
                )
	  phaseMarker[X] = 3.0f;
	else
	  phaseMarker[X] = 0.0f;
      }
    else
      phaseMarker[X] = 0.0f;
    /* >>>>>>>>> phase marking logic <<<<<<<<  */
    
  } // onsites block ends here

  //hila::out0 << "Phases marking done \n";
  
} // phaseMarking() function ends here


#define USE_PARIO
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

  
void glsol::write_moduli() {

   // real_t a = scaleFactor(t);

    double Amod = 0.0/*, AGap = 0.0*/;
    double pimod = 0.0;

    //real_t Tp[2];

    hila::set_allreduce(false);
    onsites (ALL) {
        Amod += A[X].norm();
        pimod += pi[X].norm();

    }
        
    //update_Tp(t, Tp);
    
    if (hila::myrank() == 0) {
        config.stream << t << " "
                      << Amod / lattice.volume() << " " << pimod / lattice.volume()
                      << " ";
    }

    hila::out0 << "Moduli done \n";
}// write_moduli ends at here


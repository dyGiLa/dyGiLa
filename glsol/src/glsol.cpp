#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
//#include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/fft.h"

#include "glsol.hpp"


glsol::glsol()
{
  hila::out0 << "------------------------------------------------------------" << "\n"
             << "-- dyGiLa 3-Dimensional p-Wave TDGL HPC Simulation Suites --" << "\n"
             << "------------------------------------------------------------" << std::endl;    
           
} // allocate() function ends here


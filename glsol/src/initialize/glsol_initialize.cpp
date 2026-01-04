#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/fft.h"

#include "glsol.hpp"
#include "matep.hpp"


void glsol::initialize() {
   
  real_t Tcp_mK = MP.Tcp_mK(config.Inip);
  hila::out0 <<" T_AB: "<< MP.tAB_RWS(config.Inip)*Tcp_mK << " mK. " << std::endl;
  
  switch (config.initialCondition) {
    
    case 0: { case_0(); break; }
    case 1: { case_1(); break; }
    case 2: { case_2(); break; }    
    case 3: { case_3(); break; }
    case 4: { case_4(); break; }    
    case 5: { case_5(); break; }    
    case 6: { case_6(); break; }
    case 7: { case_7(); break; } 
    case 8: { case_8(); break; } 
    case 9: { case_9(); break; } 
    
    default:
      {
       // #pragma hila ast_dump
       pi = 0.0; //set derivative matrix to zero
       deltaPi =0;
       djAaj = 0;
     
       onsites (ALL) { A[X].fill(1.0); }
    
       hila::out0 << "Field matrix set to 1 everywhere! This is default conf. " << std::endl;

       break;
      } // default block
  } // switch block ends here

} // initialize() call end here


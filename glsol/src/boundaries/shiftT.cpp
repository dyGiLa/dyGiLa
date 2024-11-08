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


const real_t glsol::shiftT(real_t Tp, real_t T0, real_t Tm, CoordinateVector p, Direction d, int dir) {

  real_t Tbc;
  CoordinateVector box;

  box={config.lx,config.ly,config.lz};

  if(p[d] == 0 and d == 0 and dir == 0){
    switch (config.bcTx0) {
    case 0: {
      Tbc = periodic_T(Tm);
      break;
    }
    } 
  }
  else if (p[d] == 0 and d == 1 and dir == 0){
    switch (config.bcTy0) {
    case 0: {
      Tbc = periodic_T(Tm);
      break;
    }
    }
  }
  else if (p[d] == 0 and d == 2 and dir == 0){
    switch (config.bcTz0) {
    case 0: {
      Tbc = periodic_T(Tm);
      break;
    }
    }
  }
  else if (p[d] == box[d]-1 and d == 0 and dir == 1){
    switch (config.bcTxN) {
    case 0: {
      Tbc = periodic_T(Tp);
      break;
    }
    }
  }
 else if (p[d] == box[d]-1 and d == 1 and dir == 1){
     switch (config.bcTyN) {
    case 0: {
      Tbc = periodic_T(Tp);
      break;
    }
    }
  }
 else if (p[d] == box[d]-1 and d == 2 and dir == 1){
   switch (config.bcTzN) {
   case 0: {
     Tbc = periodic_T(Tp);
     break;
   }
   }
 }
 else {
   if (dir == 0){
     Tbc=Tm;}
   else if (dir == 1){
     Tbc=Tp;
   }
 }
    
    
  return Tbc;

}

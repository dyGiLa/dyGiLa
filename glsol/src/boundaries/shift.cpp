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

phi_t glsol::shift(phi_t phip, phi_t phi0, phi_t phim, CoordinateVector p, Direction d, int dir) {

  phi_t phibc;
  CoordinateVector box;

  box={config.lx,config.ly,config.lz};

  if(p[d] == 0 and d == 0 and dir == 0){
    switch (config.bcx0) {
    case 0: {
      phibc = periodic(phim);
      break;
    }
    case 1: {
      phibc = AphaseBoundary();
      break;
    }
    case 2: {
      phibc = BphaseBoundary();
      break;
    }
    case 3: {
      phibc = normalBoundary();
      break;
    }
    case 4: {
      phibc = maximal();
      break;
    }
    case 5: {
      phibc = RobinBoundary(phip,0,config.robinbt_x0,0);
      break;
    }
    } 
  }
  else if (p[d] == 0 and d == 1 and dir == 0){
    switch (config.bcy0) {
    case 0: {
      phibc = periodic(phim);
      break;
    }
    case 1: {
      phibc = AphaseBoundary();
      break;
    }
    case 2: {
      phibc = BphaseBoundary();
      break;
    }
    case 3: {
      phibc = normalBoundary();
      break;
    }
    case 4: {
      phibc = maximal();
      break;
    }
      case 5: {
	phibc = RobinBoundary(phip,0,config.robinbt_y0,1);
      break;
    }
    }
  }
  else if (p[d] == 0 and d == 2 and dir == 0){
    switch (config.bcz0) {
    case 0: {
      phibc = periodic(phim);
      break;
    }
    case 1: {
      phibc = AphaseBoundary();
      break;
    }
    case 2: {
      phibc = BphaseBoundary();
      break;
    }
    case 3: {
      phibc = normalBoundary();
      break;
    }
    case 4: {
      phibc = maximal();
      break;
    }
      case 5: {
	phibc = RobinBoundary(phip,0,config.robinbt_z0,2);
      break;
    }
    }
  }
  else if (p[d] == box[d]-1 and d == 0 and dir == 1){
     switch (config.bcxN) {
    case 0: {
      phibc = periodic(phip);
      break;
    }
    case 1: {
      phibc = AphaseBoundary();
      break;
    }
    case 2: {
      phibc = BphaseBoundary();
      break;
    }
    case 3: {
      phibc = normalBoundary();
      break;
    }
    case 4: {
      phibc = maximal();
      break;
    }
      case 5: {
	phibc = RobinBoundary(phim,1,config.robinbt_xN,0);
      break;
    }
    }
  }
 else if (p[d] == box[d]-1 and d == 1 and dir == 1){
     switch (config.bcyN) {
    case 0: {
      phibc = periodic(phip);
      break;
    }
    case 1: {
      phibc = AphaseBoundary();
      break;
    }
    case 2: {
      phibc = BphaseBoundary();
      break;
    }
    case 3: {
      phibc = normalBoundary();
      break;
    }
    case 4: {
      phibc = maximal();
      break;
    }
      case 5: {
	phibc = RobinBoundary(phim,1,config.robinbt_yN,1);
      break;
    }
    }
  }
  else if (p[d] == box[d]-1 and d == 2 and dir == 1){
     switch (config.bczN) {
    case 0: {
      phibc = periodic(phip);
      break;
    }
    case 1: {
      phibc = AphaseBoundary();
      break;
    }
    case 2: {
      phibc = BphaseBoundary();
      break;
    }
    case 3: {
      phibc = normalBoundary();
      break;
    }
    case 4: {
      phibc = maximal();
      break;
    }
      case 5: {
	phibc = RobinBoundary(phim,1,config.robinbt_zN,2);
      break;
    }
    }
  }
  else {
    if (dir == 0 ){
      phibc=phim;}
    else if (dir == 1){
      phibc=phip;
    }
  }
    
    
  return phibc;

}


phi_t glsol::shift_talk(phi_t phip, phi_t phi0, phi_t phim, CoordinateVector p, Direction d, int dir) {

  phi_t phibc;
  CoordinateVector box;

  box={config.lx,config.ly,config.lz};

  if(p[d] == 0 and d == 0 and dir == 0){
    switch (config.bcx0) {
    case 0: {
      phibc = periodic(phim);
      break;
    }
    case 1: {
      phibc = AphaseBoundary();
      break;
    }
    case 2: {
      phibc = BphaseBoundary();
      break;
    }
    case 3: {
      phibc = normalBoundary();
      break;
    }
    case 4: {
      phibc = maximal();
      break;
    }
    case 5: {
      phibc = RobinBoundary(phip,0,config.robinbt_x0,0);
      break;
    }
    }
  }
  else if (p[d] == 0 and d == 1 and dir == 0){
    switch (config.bcy0) {
    case 0: {
      phibc = periodic(phim);
      hila::out0<<"periodic- \n";
      break;
    }
    case 1: {
      phibc = AphaseBoundary();
      break;
    }
    case 2: {
      phibc = BphaseBoundary();
      break;
    }
    case 3: {
      phibc = normalBoundary();
      break;
    }
    case 4: {
      phibc = maximal();
      break;
    }
      case 5: {
        phibc = RobinBoundary(phip,0,config.robinbt_y0,1);
      break;
    }
    }
  }
  else if (p[d] == 0 and d == 2 and dir == 0){
    switch (config.bcz0) {
    case 0: {
      phibc = periodic(phim);
      break;
    }
    case 1: {
      phibc = AphaseBoundary();
      break;
    }
    case 2: {
      phibc = BphaseBoundary();
      break;
    }
    case 3: {
      phibc = normalBoundary();
      break;
    }
    case 4: {
      phibc = maximal();
      break;
    }
      case 5: {
        phibc = RobinBoundary(phip,0,config.robinbt_z0,2);
      break;
    }
    }
  }
  else if (p[d] == box[d]-1 and d == 0 and dir == 1){
     switch (config.bcxN) {
    case 0: {
      phibc = periodic(phip);
      break;
    }
    case 1: {
      phibc = AphaseBoundary();
      break;
    }
    case 2: {
      phibc = BphaseBoundary();
      break;
    }
    case 3: {
      phibc = normalBoundary();
      break;
    }
    case 4: {
      phibc = maximal();
      break;
    }
      case 5: {
        phibc = RobinBoundary(phim,1,config.robinbt_xN,0);
      break;
    }
    }
  }
 else if (p[d] == box[d]-1 and d == 1 and dir == 1){
     switch (config.bcyN) {
    case 0: {
      phibc = periodic(phip);
      hila::out0<<"periodic+ \n";
      break;
    }
    case 1: {
      phibc = AphaseBoundary();
      break;
    }
    case 2: {
      phibc = BphaseBoundary();
      break;
    }
    case 3: {
      phibc = normalBoundary();
      break;
    }
    case 4: {
      phibc = maximal();
      break;
    }
      case 5: {
        phibc = RobinBoundary(phim,1,config.robinbt_yN,1);
      break;
    }
    }
  }
  else if (p[d] == box[d]-1 and d == 2 and dir == 1){
     switch (config.bczN) {
    case 0: {
      phibc = periodic(phip);
      break;
    }
    case 1: {
      phibc = AphaseBoundary();
      break;
    }
    case 2: {
      phibc = BphaseBoundary();
      break;
    }
    case 3: {
      phibc = normalBoundary();
      break;
    }
    case 4: {
      phibc = maximal();
      break;
    }
      case 5: {
        phibc = RobinBoundary(phim,1,config.robinbt_zN,2);
      break;
    }
    }
  }
  else {
    if (dir == 0 ){
      phibc=phim;}
    else if (dir == 1){
      phibc=phip;
    }
  }


  return phibc;

}

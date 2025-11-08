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


void glsol::phaseCounting() {
  
  ReductionVector<float> px_acc(/*pxacc::*/N_PMREDUCTION);
  px_acc = 0.f;
  px_acc.allreduce(true);

  // hila::set_allpxacce(false);
  onsites(ALL) {
      
      /* accumulate values */
      px_acc[/*pxacc::*/p0_acc] += (phaseMarker[X]==0) ? 1 : 0;    
      px_acc[/*pxacc::*/p1_acc] += (phaseMarker[X]==1) ? 1 : 0;
      px_acc[/*pxacc::*/p2_acc] += (phaseMarker[X]==2) ? 1 : 0;
      px_acc[/*pxacc::*/p3_acc] += (phaseMarker[X]==3) ? 1 : 0;
      px_acc[/*pxacc::*/p4_acc] += (phaseMarker[X]==4) ? 1 : 0;
      px_acc[/*pxacc::*/p5_acc] += (phaseMarker[X]==5) ? 1 : 0;
      px_acc[/*pxacc::*/p6_acc] += (phaseMarker[X]==6) ? 1 : 0;
      px_acc[/*pxacc::*/p7_acc] += (phaseMarker[X]==7) ? 1 : 0;
      px_acc[/*pxacc::*/p8_acc] += (phaseMarker[X]==8) ? 1 : 0;
      px_acc[/*pxacc::*/p9_acc] += (phaseMarker[X]==9) ? 1 : 0;

  } // onsites(ALL) block 

  std::initializer_list<int> coordsList {0,0,0};
  const CoordinateVector originpoints(coordsList);
  real_t T000 = T.get_element(originpoints);


  double vol = lattice.volume();
  float acc_vol = px_acc[/*pxacc::*/p0_acc] + px_acc[/*pxacc::*/p1_acc]
                  + px_acc[/*pxacc::*/p2_acc] + px_acc[/*pxacc::*/p3_acc]
                  + px_acc[/*pxacc::*/p4_acc] + px_acc[/*pxacc::*/p5_acc]
                  + px_acc[/*pxacc::*/p6_acc] + px_acc[/*pxacc::*/p7_acc]
                  + px_acc[/*pxacc::*/p8_acc] + px_acc[/*pxacc::*/p9_acc];

  // Volume element in unit of \xi_GL^0
  const double Velem = config.dx*config.dx*config.dx;
  // Hot Blob initial radius
  const float rb = MP.r_TcMax_blob(config.Inip, config.use_CustomerDctxi, config.Dctxi, config.Ttdb1, config.Ttdb0, config.t1) * (MP.xi0GLp(config.Inip)) * (1e6);
  
  config.streampc
         << rb << "," << t << "," << T000 << ","
	 /***************************/	 	 	 
	 << px_acc[/*pxacc::*/p0_acc]/vol << ","    
	 /***************************/	 	 	 
	 << px_acc[/*pxacc::*/p1_acc]/vol << ","
	 /***************************/	 	 
	 << px_acc[/*pxacc::*/p2_acc]/vol << ","
	 /***************************/	 
	 << px_acc[/*pxacc::*/p3_acc]/vol << "," 
	 /***************************/
	 << px_acc[/*pxacc::*/p4_acc]/vol << "," 
	 /***************************/	 
	 << px_acc[/*pxacc::*/p5_acc]/vol << "," 
	 /***************************/
	 << px_acc[/*pxacc::*/p6_acc]/vol << ","
	 /***************************/	 	 
	 << px_acc[/*pxacc::*/p7_acc]/vol << ","
	 /***************************/	 
	 << px_acc[/*pxacc::*/p8_acc]/vol << "," 
	 /***************************/
	 << px_acc[/*pxacc::*/p9_acc]/vol << ","
	 /***************************/
	 << acc_vol/vol << ","
	 /***************************/
	 << px_acc[/*pxacc::*/p0_acc]*Velem << ","    
	 /***************************/
	 << px_acc[/*pxacc::*/p1_acc]*Velem << ","
	 /***************************/	 	 
	 << px_acc[/*pxacc::*/p2_acc]*Velem << ","
	 /***************************/	 
	 << px_acc[/*pxacc::*/p3_acc]*Velem << "," 
	 /***************************/
	 << px_acc[/*pxacc::*/p4_acc]*Velem << "," 
	 /***************************/	 
	 << px_acc[/*pxacc::*/p5_acc]*Velem << "," 
	 /***************************/
	 << px_acc[/*pxacc::*/p6_acc]*Velem << ","
	 /***************************/	 	 
	 << px_acc[/*pxacc::*/p7_acc]*Velem << ","
	 /***************************/	 
	 << px_acc[/*pxacc::*/p8_acc]*Velem << "," 
	 /***************************/
	 << px_acc[/*pxacc::*/p9_acc]*Velem << ","
	 /***************************/
         << acc_vol*Velem
	 /***************************/    
         << std::endl;

  
} // phases_counter function ends here


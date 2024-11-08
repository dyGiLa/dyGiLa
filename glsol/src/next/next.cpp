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

void glsol::next() {

  static hila::timer next_timer("timestep");
  Field<phi_t> deltaPi;
  Field<Vector<3,Complex<real_t>>> djAaj;

  int bc=config.boundaryConditions;
  
  next_timer.start();

  /* --------------------------------------------------------------- */  
  /* >>>>>>>>>>  boundary condition handling block starts <<<<<<<<<< */
  /* --------------------------------------------------------------- */
  
  onsites(ALL) {

    real_t gapa = MP.gap_A_td(p, T[X]);
    real_t gapb = MP.gap_B_td(p, T[X]);

    
    A[X] += config.dt * pi[X];

    if(bc < 3)
      {
        if (bc == 1)
          {
	    if (X.coordinate(e_z) == 0 or X.coordinate(e_z) == 1)
	      {	
	       foralldir(d1)foralldir(d2){
	       if (d1==d2){
		A[X].e(d1,d2).re = 1.0;
		A[X].e(d1,d2).im = 0.0;
	       }
	       else {
		 A[X].e(d1,d2).re = 0.0;
		 A[X].e(d1,d2).im = 0.0;
	        }
	       }
	       A[X] = gapb * A[X]/sqrt(3.0);
	      }
	    else if (X.coordinate(e_z) == (config.lz - 1) or X.coordinate(e_z) == (config.lz - 2))
	        {
	         foralldir(d1)foralldir(d2){
	         if (d1==2 && d2==0){
		   A[X].e(d1,d2).re = 1.0;
		   A[X].e(d1,d2).im = 0.0;
	         }
	         else if (d1==2 && d2==1){
		   A[X].e(d1,d2).re = 0.0;  // this A-order parameter same with GL-theory note eq.46
		   A[X].e(d1,d2).im = 1.0;
	         }
	         else {
		   A[X].e(d1,d2).re = 0.0;
		   A[X].e(d1,d2).im = 0.0;
	         }
	        }
	        A[X] = gapa * A[X]/sqrt(2.0);
	        }
           }
         else if (bc == 2)
             {
	       if (
		   X.coordinate(e_x) == 0 or X.coordinate(e_x) == (config.lx - 1) or
	           X.coordinate(e_x) == 1 or X.coordinate(e_x) == (config.lx - 2) or
	           X.coordinate(e_y) == 0 or X.coordinate(e_y) == (config.ly - 1) or
	           X.coordinate(e_y) == 1 or X.coordinate(e_y) == (config.ly - 2) or
	           X.coordinate(e_z) == 0 or X.coordinate(e_z) == (config.lz - 1) or
	           X.coordinate(e_z) == 1 or X.coordinate(e_z) == (config.lz - 2)
		  )
                 {
	          A[X]=0.0;
	         }
              }

      }
    else if(bc == 3)
      {
       	if (X.coordinate(e_y) == 0 or X.coordinate(e_y) == (config.ly - 1) or
	    X.coordinate(e_y) == 1 or X.coordinate(e_y) == (config.ly - 2))
            {A[X]=0.0;}
      }
    else if(bc == 4)
      {
       	if (X.coordinate(e_y) == 0 or X.coordinate(e_y) == (config.ly - 1) or
	    X.coordinate(e_y) == 1 or X.coordinate(e_y) == (config.ly - 2) or
            X.coordinate(e_z) == 0 or X.coordinate(e_z) == (config.lz - 1) or
            X.coordinate(e_z) == 1 or X.coordinate(e_z) == (config.lz - 2))
            {A[X]=0.0;}
      }

    /*-------------- B-B domain wall BC --------------*/
    else if(bc == 5)
      {
	// if (X.coordinate(e_x) == 0 or X.coordinate(e_x) == 1)
	//       {	
	//        foralldir(d1)foralldir(d2)
	// 	 {
	//           if (d1!=d2)
	// 	  {
	// 	   A[X].e(d1,d2).re = 0.0;
	// 	   A[X].e(d1,d2).im = 0.0;
	//           }
	//           else
	// 	  {
	// 	   A[X].e(0,0).re = config.BLeft_11;
	// 	   A[X].e(0,0).im = 0.0;
	// 	   A[X].e(1,1).re = config.BLeft_22;
	// 	   A[X].e(1,1).im = 0.0;
	// 	   A[X].e(2,2).re = config.BLeft_33;
	// 	   A[X].e(2,2).im = 0.0;		   
	//           }
 	//           A[X] = gapb * A[X]/sqrt(3.0); 
	//          } // foralldir ends here
	//       }
	//     else if (X.coordinate(e_x) == (config.lx - 1) or X.coordinate(e_x) == (config.lx - 2))
	//         {
	//          foralldir(d1)foralldir(d2)
	// 	 {
	//           if (d1!=d2)
	// 	  {
	// 	   A[X].e(d1,d2).re = 0.0;
	// 	   A[X].e(d1,d2).im = 0.0;
	//           }
	//           else
	// 	  {
	//            A[X].e(0,0).re = config.BRight_11;
	// 	   A[X].e(0,0).im = 0.0;
	// 	   A[X].e(1,1).re = config.BRight_22;
	// 	   A[X].e(1,1).im = 0.0;
	// 	   A[X].e(2,2).re = config.BRight_33;
	// 	   A[X].e(2,2).im = 0.0;		   
	//           }
  	//           A[X] = gapb * A[X]/sqrt(3.0);
	//          }
	//         } // foralldir ends here
      }
      /*---------- B-B domain wall BC ends here ----------*/

      /*--------- B-phase phaseVortices BC   -------------*/
    else if(bc==6)
      {
        // if (
	//     (X.coordinate(e_x) == 0 || X.coordinate(e_x) == (config.lx - 1))
	//     || (X.coordinate(e_y) == 0 || X.coordinate(e_y) == (config.ly - 1))
	//    )
	//   {
	//     real_t mod       = sqrt((X.coordinate(e_x)-(config.lx)/2.) * (X.coordinate(e_x)-(config.lx)/2.)
	// 		            + (X.coordinate(e_y)-(config.ly)/2.) * (X.coordinate(e_y)-(config.ly)/2.));
	//     Complex<real_t> phaseExp; // exp(i \phi)
	//     phaseExp.re = (X.coordinate(e_x)-(config.lx)/2.)/mod;
	//     phaseExp.im = (X.coordinate(e_y)-(config.ly)/2.)/mod;
	    
	//     foralldir(i) foralldir(al)
	//       {
        //         if (i != al)
	// 	  A[X].e(i,al) = 0.;
	// 	else
	// 	  {
	// 	    /*A[X].e(i,al).re = (gapb/sqrt(3.)) * ((X.coordinate(e_x)-(config.lx)/2.)/mod);
	// 	      A[X].e(i,al).im = (gapb/sqrt(3.)) * ((X.coordinate(e_y)-(config.ly)/2.)/mod);*/
	// 	    A[X].e(i,al).re = gapb/sqrt(3.);
	// 	    A[X].e(i,al).im = 0.;		    
	// 	  }
	//       } // SO(3) R of A
	    
	//     for (unsigned int n = 0; n<config.Wn; ++n) {A[X] = A[X] * phaseExp;} // R e^{Wn\phi}
	//   }
	  
      } // bc = 6 block, phase vortices ends here
    /*------------    phaseVortices BC ends ---------------*/
    
  } // onsites(ALL) block ends here

  /* --------------------------------------------------------------- */  
  /* >>>>>>>    boundary condition handling block ends here   <<<<<< */
  /* --------------------------------------------------------------- */  

  onsites (ALL) {

    real_t beta[6];
    point_params(T[X], p,beta);
      
    auto AxAt = A[X]*A[X].transpose();
    auto AxAd = A[X]*A[X].dagger();
    
    deltaPi[X] = - beta[0]*A[X]
      - 2.0*beta[1]*A[X].conj()*AxAt.trace() 
      - 2.0*beta[2]*A[X]*AxAd.trace()
      - 2.0*beta[3]*AxAt*A[X].conj()
      - 2.0*beta[4]*AxAd*A[X] 
      - 2.0*beta[5]*A[X].conj()*A[X].transpose()*A[X];

  } // bulk energy contribution

  onsites(ALL) {
    djAaj[X] = 0;
    foralldir(j) {
      djAaj[X] += A[X + j].column(j) - A[X - j].column(j);
    }
  } // DjA_aj = D0A_al0 + D1A_al1 + D2A_al2 ???
  
  onsites(ALL) {
    phi_t mat;
    foralldir(d) {
      auto col = djAaj[X+d] - djAaj[X-d];
      for (int i=0; i<NDIM; i++) mat.e(i,d) = col[i];
    }

    deltaPi[X] += (1.0/(2.0*(config.dx*config.dx)))*mat;
  }
  
  onsites(ALL) {

    deltaPi[X] +=  (1.0 / (config.dx * config.dx)) * (A[X + e_x] + A[X - e_x]
						      + A[X + e_y] + A[X - e_y]
						      + A[X + e_z] + A[X - e_z]
						      - 6.0 * A[X]);
  } // DjDjAali term summation

  
  if (t < config.tdif)
    {
      pi[ALL] = deltaPi[X]/(config.difFac);
      t += config.dt/config.difFac;
    }
  else if (t < config.tdis && config.gamma.squarenorm() > 0 )
    {
      //hila::out0 << "config.gamma is " << config.gamma << "\n" << std::endl;
      pi[ALL] = pi[X] + (deltaPi[X] - 2.0 * config.gamma * pi[X])*config.dt; //Complex<real_t> C(a, b) = r + I *
      t += config.dt;
    }
  else
    {
      pi[ALL] = pi[X] + deltaPi[X]*config.dt;
      t += config.dt;
    }

  next_timer.stop();

} // next() function ends here


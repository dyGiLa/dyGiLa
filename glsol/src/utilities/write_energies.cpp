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


void glsol::write_energies() {

  // Complex<double> suma(0),sumb2(0),sumb3(0),sumb4(0),sumb5(0);
  // Complex<double> suma_we(0),sumb2_we(0),sumb3_we(0),sumb4_we(0),sumb5_we(0);
  // Complex<double> sumk1(0),sumk2(0),sumk3(0);
  // Complex<double> sumk1_we(0),sumk2_we(0),sumk3_we(0);
  // real_t sumb1 = 0;
  // Complex<double> sumb1_we = 0;
  // Complex<double> sumkin(0);
  // Complex<double> sumkin_we(0);
  // real_t Tp[2];

  // update_Tp(t, Tp);

  //real_t ebfe=fmin(MP.f_A_td(Tp[1], Tp[0]),MP.f_B_td(Tp[1], Tp[0]));

  ReductionVector<Complex<double>> red(N_REDUCTION);
  red = 0;
  red.allreduce(false);

  // hila::set_allreduce(false);
  onsites(ALL) {

      Matep MPonsites;
      Complex<double> a(0),b2(0),b3(0),b4(0),b5(0);
      Complex<double> kin(0);
      Complex<double> k1(0), k2(0), k3(0);
      Complex<double> bfe(0);
      double b1 = 0;

      real_t ebfe = MPonsites.f_A_td(p[X], T[X]);// fmin(MP.f_A_td(p[X], T[X]),MP.f_B_td(p[X], T[X]));
      
      // local array of alpha, beta_i, alpha = beta[0]
      real_t beta[6];
      beta[0] = MPonsites.alpha_td(p[X], T[X]);
      beta[1] = MPonsites.beta1_td(p[X], T[X]);
      beta[2] = MPonsites.beta2_td(p[X], T[X]);
      beta[3] = MPonsites.beta3_td(p[X], T[X]);
      beta[4] = MPonsites.beta4_td(p[X], T[X]);
      beta[5] = MPonsites.beta5_td(p[X], T[X]);
      // point_params(T[X], p[X], beta);
      
      a = beta[0] * (A[X]*A[X].dagger()).trace();

      b1 = beta[1] * ((A[X]*A[X].transpose()).trace()).squarenorm();

      b2 = beta[2] * ((A[X]*A[X].dagger()).trace()*(A[X]*A[X].dagger()).trace());

      b3 = beta[3] * ((A[X]*A[X].transpose()*A[X].conj()*A[X].dagger()).trace());

      b4 = beta[4] * ((A[X]*A[X].dagger()*A[X]*A[X].dagger()).trace());

      b5 = beta[5] * ((A[X]*A[X].dagger()*A[X].conj()*A[X].transpose()).trace());

      bfe = a + b1 + b2 + b3 + b4 + b5 - ebfe;

      kin = (pi[X]*pi[X].dagger()).trace();
      
      foralldir(j) foralldir (k) foralldir(al){
          k1 += squarenorm(A[X + k].e(al,j) - A[X - k].e(al,j)) / (4.0*sqr(config.dx));
          k2 += (A[X + j].e(al,j) - A[X - j].e(al,j)) * (A[X + k].e(al,k) - A[X - k].e(al,k)).conj() / (4.0*sqr(config.dx));
          k3 += (A[X + k].e(al,j) - A[X - k].e(al,j)) * (A[X + j].e(al,k) - A[X - j].e(al,k)).conj() / (4.0*sqr(config.dx));
      } // block of gridient energy terms 

      red[i_sumkin] += kin;
      red[i_sumkin_we] += bfe * kin;

      red[i_sumk1] += k1;
      red[i_sumk1_we] += bfe * k1;

      red[i_sumk2] += k2;
      red[i_sumk2_we] += bfe *	k2;

      red[i_sumk3] += k3;
      red[i_sumk3_we] += bfe *	k3;
      
      red[i_suma] += a;
      red[i_suma_we] += bfe * a;

      // sumb1 += b1;
      red[i_sumb1] += b1;
      red[i_sumb1_we] += bfe * b1;

      red[i_sumb2] += b2;
      red[i_sumb2_we] += bfe * b2;

      red[i_sumb3] += b3;
      red[i_sumb3_we] += bfe * b3;

      red[i_sumb4] += b4;
      red[i_sumb4_we] += bfe * b4;

      red[i_sumb5] += b5;
      red[i_sumb5_we] += bfe * b5;
  } // onsites(ALL) block 

  if (hila::myrank() == 0) {
       double vol = lattice.volume();
       config.stream 
	 << red[i_sumkin].re / vol << " " << red[i_sumkin].im / vol << " "
         << red[i_sumkin_we].re / vol << " " << red[i_sumkin_we].re / vol << " "
	 /***************************/
         << red[i_sumk1].re / vol << " " << red[i_sumk1].im / vol << " "
         << red[i_sumk1_we].re / vol << " " << red[i_sumk1_we].im / vol << " "
	 /***************************/	 
         << red[i_sumk2].re / vol << " " << red[i_sumk2].im / vol << " "
         << red[i_sumk2_we].re / vol << " " << red[i_sumk2_we].im / vol << " "
	 /***************************/	 
         << red[i_sumk3].re / vol << " " << red[i_sumk3].im / vol << " "
         << red[i_sumk3_we].re / vol << " " << red[i_sumk3_we].im / vol << " "
	 /***************************/	 
         << red[i_suma].re / vol << " " << red[i_suma].im / vol << " "
         << red[i_suma_we].re / vol << " " << red[i_suma_we].im / vol << " "
	 /***************************/	 
         << red[i_sumb1].re / vol << " " << red[i_sumb1].im / vol << " "
         << red[i_sumb1_we].re / vol << " " << red[i_sumb1_we].im / vol << " "
	 /***************************/	 
         << red[i_sumb2].re / vol << " " << red[i_sumb2].im / vol << " "
         << red[i_sumb2_we].re / vol << " " << red[i_sumb2_we].im / vol << " "
	 /***************************/	 
         << red[i_sumb3].re / vol << " " << red[i_sumb3].im / vol << " " 
         << red[i_sumb3_we].re / vol << " " << red[i_sumb3_we].im / vol << " "
	 /***************************/	 
         << red[i_sumb4].re / vol << " " << red[i_sumb4].im / vol << " "
         << red[i_sumb4_we].re / vol << " " << red[i_sumb4_we].im / vol << " "
	 /***************************/	 
         << red[i_sumb5].re / vol << " " << red[i_sumb5].im / vol << " "
         << red[i_sumb5_we].re / vol << " " << red[i_sumb5_we].im / vol << " "
	 /***************************/	 
         << std::endl;

  } // config.stream output block

  hila::out0 << "Energy done \n";
  
} // write_energies function ends here


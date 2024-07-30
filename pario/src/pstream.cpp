#define _USE_MATH_DEFINES
#define USE_ASCENT 
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
#include "pario.hpp"

//#if defined USE_ASCENT
#include "ascent.hpp"
#include "conduit_blueprint.hpp"
//#endif

void parIO::pstream(glsol &sol) {

    /*-------------------    sqrt(Tr[A.A^+) ----------------------*/
    gapA[ALL] = real(sqrt((sol.A[X]*(sol.A[X].dagger())).trace()));  


    // /*--------------------     feDensity      --------------------*/
    // real_t ebfe=fmin(matep.f_A_td(config.p, config.T), matep.f_B_td(config.p, config.T));
    // feDensity[ALL] = 0;
    // onsites(ALL) {
    //   Complex<double> a(0),b2(0),b3(0),b4(0),b5(0);
    //   //Complex<double> kin(0);
    //   Complex<double> k1(0), k2(0), k3(0);
    //   Complex<double> bfe(0);
    //   double b1 = 0;
      
    //   a = config.alpha * (A[X]*A[X].dagger()).trace();
    //   b1 = config.beta1 * ((A[X]*A[X].transpose()).trace()).squarenorm();
    //   b2 = config.beta2 * ((A[X]*A[X].dagger()).trace()*(A[X]*A[X].dagger()).trace());
    //   b3 = config.beta3 * ((A[X]*A[X].transpose()*A[X].conj()*A[X].dagger()).trace());
    //   b4 = config.beta4 * ((A[X]*A[X].dagger()*A[X]*A[X].dagger()).trace());
    //   b5 = config.beta5 * ((A[X]*A[X].dagger()*A[X].conj()*A[X].transpose()).trace());

    //   bfe = a + b1 + b2 + b3 + b4 + b5 - ebfe;
    //   //kin = (pi[X]*pi[X].dagger()).trace();
      
    //   foralldir(j) foralldir (k) foralldir(al){
    // 	k1 += (A[X + k].column(j) - A[X - k].column(j)).e(al)
    // 	      * (A[X + k].conj().column(j) - A[X - k].conj().column(j)).e(al)/(4.0*config.dx*config.dx);
    // 	k2 += (A[X + j].column(j) - A[X - j].column(j)).e(al)
    // 	      * (A[X + k].conj().column(k) - A[X - k].conj().column(k)).e(al)/(4.0*config.dx*config.dx);
    // 	k3 += (A[X + k].column(j) - A[X - k].column(j)).e(al)
    // 	      * (A[X + j].conj().column(k) - A[X - j].conj().column(k)).e(al)/(4.0*config.dx*config.dx);
    //   }
    //   // question here : how about imagnary part of k1 + k2 + k3 +bfe
    //   feDensity[X] = real(k1 + k2 + k3 + bfe);
    // } //onsite(All) end here
    
    // /*----------------     A matrix elements ---------------------*/
    // if (config.A_matrix_output == 1){
    //  u11[ALL] = A[X].e(0,0).re; v11[ALL] = A[X].e(0,0).im;
    //  u12[ALL] = A[X].e(0,1).re; v12[ALL] = A[X].e(0,1).im;
    //  u13[ALL] = A[X].e(0,2).re; v13[ALL] = A[X].e(0,2).im;
    //  u21[ALL] = A[X].e(1,0).re; v21[ALL] = A[X].e(1,0).im;
    //  u22[ALL] = A[X].e(1,1).re; v22[ALL] = A[X].e(1,1).im;
    //  u23[ALL] = A[X].e(1,2).re; v23[ALL] = A[X].e(1,2).im;
    //  u31[ALL] = A[X].e(2,0).re; v31[ALL] = A[X].e(2,0).im;
    //  u32[ALL] = A[X].e(2,1).re; v32[ALL] = A[X].e(2,1).im;
    //  u33[ALL] = A[X].e(2,2).re; v33[ALL] = A[X].e(2,2).im;

    //  u11.copy_local_data_with_halo(u11Ordered); v11.copy_local_data_with_halo(v11Ordered);
    //  u12.copy_local_data_with_halo(u12Ordered); v12.copy_local_data_with_halo(v12Ordered);
    //  u13.copy_local_data_with_halo(u13Ordered); v13.copy_local_data_with_halo(v13Ordered);
    //  u21.copy_local_data_with_halo(u21Ordered); v21.copy_local_data_with_halo(v21Ordered);
    //  u22.copy_local_data_with_halo(u22Ordered); v22.copy_local_data_with_halo(v22Ordered);
    //  u23.copy_local_data_with_halo(u23Ordered); v23.copy_local_data_with_halo(v23Ordered);
    //  u31.copy_local_data_with_halo(u31Ordered); v31.copy_local_data_with_halo(v31Ordered);
    //  u32.copy_local_data_with_halo(u32Ordered); v32.copy_local_data_with_halo(v32Ordered);
    //  u33.copy_local_data_with_halo(u33Ordered); v33.copy_local_data_with_halo(v33Ordered);        
    // }

    // /*----------------------     trace A    ----------------------*/
    // if (config.hdf5_trA_output == 1){
    //   trA_re[ALL] = (A[X].trace()).real();
    //   trA_im[ALL] = (A[X].trace()).imag();      
    //   trA_re.copy_local_data_with_halo(trA_reOrdered);
    //   trA_im.copy_local_data_with_halo(trA_imOrdered);      
    // }

    // /*----------------------  eigen values of A ------------------*/
    // if (config.hdf5_eigvA_output == 1){
    //   Field<Vector<3,double>> eval;
    //   Field<Matrix<3,3,Complex<double>>> evec;

    //   onsites(ALL){
    // 	A[X].eigen_jacobi(eval[X],evec[X]/*,hila::sort::ascending*/);
    //   }

    //   eigAv1[ALL] = eval[X].e(0); 
    //   eigAv2[ALL] = eval[X].e(1); 
    //   eigAv3[ALL] = eval[X].e(2); 

    //   eigAv1.copy_local_data_with_halo(eigAv1Ordered);
    //   eigAv2.copy_local_data_with_halo(eigAv2Ordered);
    //   eigAv3.copy_local_data_with_halo(eigAv3Ordered);
    // }

    // /*------------------ mass current components ------------------*/
    // if (config.hdf5_mass_current_output == 1){
    //   Field<Vector<3,double>> jmX;
    //   Field<Complex<real_t>>  phaseExp;

    //   /*onsites(ALL){
    //     foralldir(i) foralldir(j) foralldir(al){
    // 	  jmX[X].e(i) = ((A[X].conj().column(j)).e(al) * (A[X+i].column(j) - A[X-i].column(j)).e(al)/(2.*config.dx)
    // 	                + (A[X].conj().column(j)).e(al) * (A[X+j].column(i) - A[X-j].column(i)).e(al)/(2.*config.dx)
    // 			+ (A[X].conj().column(i)).e(al) * (A[X+j].column(j) - A[X-j].column(j)).e(al)/(2.*config.dx)).imag();

    //     } // foralldir() calls end here

    // 	} //onesite(ALL) call ends here*/

    //   onsites(ALL) {
    // 	jmX[X] = 0;
    // 	foralldir(i) foralldir(j) foralldir(al) {
    // 	  jmX[X].e(i) += (A[X].e(al,j).conj() *(A[X+i].e(al,j) - A[X-i].e(al,j))
    // 			  + A[X].e(al,j).conj() * (A[X+j].e(al,i) - A[X-j].e(al,i))
    // 			  + A[X].e(al,i).conj() * (A[X+j].e(al,j) - A[X-j].e(al,j))).imag();
    // 	} // foralldir end here, outermost foralldir slowest, inner run earier
    // 	jmX[X] /= 2*config.dx;
    //   } // onsites(ALL) end here

    //   jm1[ALL] = jmX[X].e(0);
    //   jm2[ALL] = jmX[X].e(1);
    //   jm3[ALL] = jmX[X].e(2);

    //   /* >>>>>>>> Modulus, phase angle, phaseExp   <<<<<<< */
    //   phaseExp[ALL]        = ((A[X].transpose()) * A[X]).trace();
    //   phaseExp2Re[ALL]    = phaseExp[X].real();
    //   phaseExp2Im[ALL]    = phaseExp[X].imag();      
      
    //   phaseExpModulus[ALL] = phaseExp[X].abs();
    //   //phaseExpAngle[ALL]   = phaseExp[X].arg()/2.;
    //   phaseExpAngle[ALL]   = std::atan2(phaseExp[X].imag(), phaseExp[X].real())/2.;
    //   /* >>>>>>>> phase angle and modulus end here   <<<<<<*/

    //   jm1.copy_local_data_with_halo(jm1Ordered);
    //   jm2.copy_local_data_with_halo(jm2Ordered);
    //   jm3.copy_local_data_with_halo(jm3Ordered);

    //   phaseExpModulus.copy_local_data_with_halo(phaseExpModulusO);
    //   phaseExpAngle.copy_local_data_with_halo(phaseExpAngleO);
    //   phaseExp2Re.copy_local_data_with_halo(phaseExp2ReO);
    //   phaseExp2Im.copy_local_data_with_halo(phaseExp2ImO);            
    // }

    // /*------------------ spin current components ------------------*/
    // if (config.hdf5_spin_current_output == 1){
    //   Field<Matrix<3,3,double>> jsX; // column is alpha for spin, row is i for spatial

    //   onsites(ALL) {
    // 	jsX[X] = 0;
    // 	foralldir(al) foralldir(i) foralldir(be) foralldir(ga) foralldir(j) {
    //       jsX[X].e(i,al) += -matep.epsilon(al,be,ga)
    // 	                     * (A[X].e(be,i).conj() * (A[X+j].e(ga,j) - A[X-j].e(ga,j))
    // 			        + A[X].e(be,j).conj() * (A[X+i].e(ga,j) - A[X-i].e(ga,j))
    // 			        + A[X].e(be,j).conj() * (A[X+j].e(ga,i) - A[X-j].e(ga,i))).real();
	  
    // 	} // foralldir end here, outermost foralldir slowest, inner run earier
    // 	jsX[X] /= 2*config.dx;
    //   } // onsites(ALL) end here

    //   js11[ALL] = jsX[X].e(0,0);
    //   js21[ALL] = jsX[X].e(1,0);
    //   js31[ALL] = jsX[X].e(2,0);
    //   js12[ALL] = jsX[X].e(0,1);
    //   js22[ALL] = jsX[X].e(1,1);
    //   js32[ALL] = jsX[X].e(2,1);
    //   js13[ALL] = jsX[X].e(0,2);
    //   js23[ALL] = jsX[X].e(1,2);
    //   js33[ALL] = jsX[X].e(2,2);
      
    //   js11.copy_local_data_with_halo(js11O);
    //   js21.copy_local_data_with_halo(js21O);
    //   js31.copy_local_data_with_halo(js31O);
      
    //   js12.copy_local_data_with_halo(js12O);
    //   js22.copy_local_data_with_halo(js22O);
    //   js32.copy_local_data_with_halo(js32O);
      
    //   js13.copy_local_data_with_halo(js13O);
    //   js23.copy_local_data_with_halo(js23O);
    //   js33.copy_local_data_with_halo(js33O);

    // }
    
    /*------------------------------------------------------------*/
        
    gapA.copy_local_data_with_halo(gapAOrdered);
    //feDensity.copy_local_data_with_halo(feDensityOrdered);
    
    pio.execute(actions);

} // pstream() end here


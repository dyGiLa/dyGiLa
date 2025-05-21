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
#include "pario.hpp"

#include "ascent.hpp"
#include "conduit_blueprint.hpp"

void parIO::pstream(glsol &sol) {

    /*-------------------    sqrt(Tr[A.A^+) ----------------------*/
    gapA[ALL] = real(sqrt((sol.A[X]*(sol.A[X].dagger())).trace()));  
    gapA.copy_local_data_with_halo(gapAContainer);

    /*--------------------     feDensity      --------------------*/
    real_t ebfe=fmin(matep.f_A_td(sol.config.Inip, sol.config.IniT), matep.f_B_td(sol.config.Inip, sol.config.IniT));
    feDensity[ALL] = 0;
    onsites(ALL) {
      Complex<double> a(0),b2(0),b3(0),b4(0),b5(0);
      //Complex<double> kin(0);
      Complex<double> k1(0), k2(0), k3(0);
      Complex<double> bfe(0);
      double b1 = 0;
      
      a = sol.config.alpha * (sol.A[X]*sol.A[X].dagger()).trace();
      b1 = sol.config.beta1 * ((sol.A[X]*sol.A[X].transpose()).trace()).squarenorm();
      b2 = sol.config.beta2 * ((sol.A[X]*sol.A[X].dagger()).trace()*(sol.A[X]*sol.A[X].dagger()).trace());
      b3 = sol.config.beta3 * ((sol.A[X]*sol.A[X].transpose()*sol.A[X].conj()*sol.A[X].dagger()).trace());
      b4 = sol.config.beta4 * ((sol.A[X]*sol.A[X].dagger()*sol.A[X]*sol.A[X].dagger()).trace());
      b5 = sol.config.beta5 * ((sol.A[X]*sol.A[X].dagger()*sol.A[X].conj()*sol.A[X].transpose()).trace());

      bfe = a + b1 + b2 + b3 + b4 + b5 - ebfe;
      //kin = (pi[X]*pi[X].dagger()).trace();
      
      foralldir(j) foralldir (k) foralldir(al){
	k1 += (sol.A[X + k].column(j) - sol.A[X - k].column(j)).e(al)
	      * (sol.A[X + k].conj().column(j) - sol.A[X - k].conj().column(j)).e(al)/(4.0*sol.config.dx*sol.config.dx);
	k2 += (sol.A[X + j].column(j) - sol.A[X - j].column(j)).e(al)
	      * (sol.A[X + k].conj().column(k) - sol.A[X - k].conj().column(k)).e(al)/(4.0*sol.config.dx*sol.config.dx);
	k3 += (sol.A[X + k].column(j) - sol.A[X - k].column(j)).e(al)
	      * (sol.A[X + j].conj().column(k) - sol.A[X - j].conj().column(k)).e(al)/(4.0*sol.config.dx*sol.config.dx);
      }
      // question here : how about imagnary part of k1 + k2 + k3 +bfe
      feDensity[X] = real(k1 + k2 + k3 + bfe);
    } //onsite(All) end here

    feDensity.copy_local_data_with_halo(feDensityContainer);

    /*-------------------- Temperature field --------------------*/
    sol.T.copy_local_data_with_halo(Temperature);

    
    /*-------------------- phaseMarker field --------------------*/
    sol.phaseMarker.copy_local_data_with_halo(phaseMarker);    

    
    /*------------------ mass current components ------------------*/
    if (sol.config.hdf5_mass_current_output == 1){
      Field<Vector<3,double>> jmX;
      Field<Complex<real_t>>  phaseExp;

      /*onsites(ALL){
        foralldir(i) foralldir(j) foralldir(al){
	  jmX[X].e(i) = ((A[X].conj().column(j)).e(al) * (A[X+i].column(j) - A[X-i].column(j)).e(al)/(2.*config.dx)
	                + (A[X].conj().column(j)).e(al) * (A[X+j].column(i) - A[X-j].column(i)).e(al)/(2.*config.dx)
			+ (A[X].conj().column(i)).e(al) * (A[X+j].column(j) - A[X-j].column(j)).e(al)/(2.*config.dx)).imag();

        } // foralldir() calls end here

	} //onesite(ALL) call ends here*/

      onsites(ALL) {
	jmX[X] = 0;
	foralldir(i) foralldir(j) foralldir(al) {
	  jmX[X].e(i) += (sol.A[X].e(al,j).conj() *(sol.A[X+i].e(al,j) - sol.A[X-i].e(al,j))
			  + sol.A[X].e(al,j).conj() * (sol.A[X+j].e(al,i) - sol.A[X-j].e(al,i))
			  + sol.A[X].e(al,i).conj() * (sol.A[X+j].e(al,j) - sol.A[X-j].e(al,j))).imag();
	} // foralldir end here, outermost foralldir slowest, inner run earier
	jmX[X] /= 2*sol.config.dx;
      } // onsites(ALL) end here

      jm1[ALL] = jmX[X].e(0);
      jm2[ALL] = jmX[X].e(1);
      jm3[ALL] = jmX[X].e(2);

      /* >>>>>>>> Modulus, phase angle, phaseExp   <<<<<<< */
      phaseExp[ALL]        = ((sol.A[X].transpose()) * sol.A[X]).trace();
      phaseExp2Re[ALL]    = phaseExp[X].real();
      phaseExp2Im[ALL]    = phaseExp[X].imag();      
      
      phaseExpModulus[ALL] = phaseExp[X].abs();
      //phaseExpAngle[ALL]   = phaseExp[X].arg()/2.;
      phaseExpAngle[ALL]   = std::atan2(phaseExp[X].imag(), phaseExp[X].real())/2.;
      /* >>>>>>>> phase angle and modulus end here   <<<<<<*/

      jm1.copy_local_data_with_halo(jm1Container);
      jm2.copy_local_data_with_halo(jm2Container);
      jm3.copy_local_data_with_halo(jm3Container);

      phaseExpModulus.copy_local_data_with_halo(phaseExpModulusContanier);
      phaseExpAngle.copy_local_data_with_halo(phaseExpAngleContanier);
      phaseExp2Re.copy_local_data_with_halo(phaseExp2ReContanier);
      phaseExp2Im.copy_local_data_with_halo(phaseExp2ImContanier);            
    }

    // /*------------------ spin current components ------------------*/
    if (sol.config.hdf5_spin_current_output == 1){
      Field<Matrix<3,3,double>> jsX; // column is alpha for spin, row is i for spatial

      onsites(ALL) {
	matep::Matep MPonsites;
	jsX[X] = 0;
	foralldir(al) foralldir(i) foralldir(be) foralldir(ga) foralldir(j) {
          jsX[X].e(i,al) += -MPonsites.epsilon(al,be,ga)
	                     * (sol.A[X].e(be,i).conj() * (sol.A[X+j].e(ga,j) - sol.A[X-j].e(ga,j))
			        + sol.A[X].e(be,j).conj() * (sol.A[X+i].e(ga,j) - sol.A[X-i].e(ga,j))
			        + sol.A[X].e(be,j).conj() * (sol.A[X+j].e(ga,i) - sol.A[X-j].e(ga,i))).real();
	  
	} // foralldir end here, outermost foralldir slowest, inner run earier
	jsX[X] /= 2*sol.config.dx;
      } // onsites(ALL) end here

      js11[ALL] = jsX[X].e(0,0);
      js21[ALL] = jsX[X].e(1,0);
      js31[ALL] = jsX[X].e(2,0);
      js12[ALL] = jsX[X].e(0,1);
      js22[ALL] = jsX[X].e(1,1);
      js32[ALL] = jsX[X].e(2,1);
      js13[ALL] = jsX[X].e(0,2);
      js23[ALL] = jsX[X].e(1,2);
      js33[ALL] = jsX[X].e(2,2);
      
      js11.copy_local_data_with_halo(js11Contanier);
      js21.copy_local_data_with_halo(js21Contanier);
      js31.copy_local_data_with_halo(js31Contanier);
      
      js12.copy_local_data_with_halo(js12Contanier);
      js22.copy_local_data_with_halo(js22Contanier);
      js32.copy_local_data_with_halo(js32Contanier);
      
      js13.copy_local_data_with_halo(js13Contanier);
      js23.copy_local_data_with_halo(js23Contanier);
      js33.copy_local_data_with_halo(js33Contanier);

    }

    /*----------------     A matrix elements ---------------------*/
    if (sol.config.hdf5_A_matrix_output == 1){
     u11[ALL] = sol.A[X].e(0,0).re; v11[ALL] = sol.A[X].e(0,0).im;
     u12[ALL] = sol.A[X].e(0,1).re; v12[ALL] = sol.A[X].e(0,1).im;
     u13[ALL] = sol.A[X].e(0,2).re; v13[ALL] = sol.A[X].e(0,2).im;
     u21[ALL] = sol.A[X].e(1,0).re; v21[ALL] = sol.A[X].e(1,0).im;
     u22[ALL] = sol.A[X].e(1,1).re; v22[ALL] = sol.A[X].e(1,1).im;
     u23[ALL] = sol.A[X].e(1,2).re; v23[ALL] = sol.A[X].e(1,2).im;
     u31[ALL] = sol.A[X].e(2,0).re; v31[ALL] = sol.A[X].e(2,0).im;
     u32[ALL] = sol.A[X].e(2,1).re; v32[ALL] = sol.A[X].e(2,1).im;
     u33[ALL] = sol.A[X].e(2,2).re; v33[ALL] = sol.A[X].e(2,2).im;

     u11.copy_local_data_with_halo(u11Container); v11.copy_local_data_with_halo(v11Container);
     u12.copy_local_data_with_halo(u12Container); v12.copy_local_data_with_halo(v12Container);
     u13.copy_local_data_with_halo(u13Container); v13.copy_local_data_with_halo(v13Container);
     u21.copy_local_data_with_halo(u21Container); v21.copy_local_data_with_halo(v21Container);
     u22.copy_local_data_with_halo(u22Container); v22.copy_local_data_with_halo(v22Container);
     u23.copy_local_data_with_halo(u23Container); v23.copy_local_data_with_halo(v23Container);
     u31.copy_local_data_with_halo(u31Container); v31.copy_local_data_with_halo(v31Container);
     u32.copy_local_data_with_halo(u32Container); v32.copy_local_data_with_halo(v32Container);
     u33.copy_local_data_with_halo(u33Container); v33.copy_local_data_with_halo(v33Container);        
    }
    
    /*------------------------------------------------------------*/
         
    pio.execute(actions);

} // pstream() end here


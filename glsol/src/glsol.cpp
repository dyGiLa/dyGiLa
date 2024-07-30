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

// /*-----------------------------------------------------------------------*/
// /***      include Ascent & Canduit for in situ rank rendering        *****/
// /*-----------------------------------------------------------------------*/
// #if defined USE_ASCENT

// #include "ascent.hpp"
// #include "conduit_blueprint.hpp"

// #endif
// /*-----------------------------------------------------------------------*/
// /***               include Ascent & Canduit end here                 *****/
// /*-----------------------------------------------------------------------*/

const std::string glsol::allocate(const std::string &fname, int argc, char **argv)
{
    Matep MP;
    hila::initialize(argc, argv);
    hila::input parameters(fname);
    
    config.lx = parameters.get("Nx");
    config.ly = parameters.get("Ny");
    config.lz = parameters.get("Nz");
    config.dx = parameters.get("dx");
    config.dtdxRatio = parameters.get("dtdxRatio");
    config.tStart = parameters.get("tStart");
    config.tEnd = parameters.get("tEnd");
    config.tdif = parameters.get("tdif");
    config.difFac = parameters.get("difFac");
    config.tdis = parameters.get("tdis");

    /*----   gamma as a complex number   -----*/
    //config.gamma = parameters.get("gamma");
    std::vector<real_t> tmp1 = parameters.get("gamma1");
    std::vector<real_t> tmp2 = parameters.get("gamma2");    

    if (tmp1.size() < 1 || tmp1.size() >2)
      {
	hila::out0 << "error: gamma1 must be initialized by at least one real or imagnary parts" << std::endl;
	hila::error("\n");
      }
    else if (tmp1.size() == 1)
      {
	hila::out0 << " tmp1.size() = " << tmp1.size() << "\n";
        config.gamma1.real() = tmp1[0];
        config.gamma1.imag() = 0.;
      }
    else  // tmp.size() == 2
      {
        hila::out0 << " tmp1.size() = " << tmp1.size() << "\n";
        config.gamma1.real() = tmp1[0];
        config.gamma1.imag() = tmp1[1];
      }

    if (tmp2.size() < 1 || tmp2.size() >2)
      {
	hila::out0 << "error: gamma must be initialized by at least one real or imagnary parts" << std::endl;
	hila::error("\n");
      }
    else if (tmp2.size() == 1)
      {
	hila::out0 << " tmp2.size() = " << tmp2.size() << "\n";
        config.gamma2.real() = tmp2[0];
        config.gamma2.imag() = 0.;
      }
    else  // tmp.size() == 2
      {
        hila::out0 << " tmp2.size() = " << tmp2.size() << "\n";
        config.gamma2.real() = tmp2[0];
        config.gamma2.imag() = tmp2[1];
      }
    
    /*---- gamma as a complex number ends-----*/
    config.gammaoffc = parameters.get("gammaoffc");

    config.initialCondition = parameters.get_item("initialCondition",{"gaussrand"             //0
								      ,"kgaussrand"           //1
								      ,"normal_phase_real1"   //2
								      ,"normal_phase_real2"   //3
								      ,"normal_phase_complex" //4
								      ,"Bphase"               //5
								      ,"Aphase_partial1"      //6
								      ,"Aphase_full"          //7
                                                                      ,"BinA"});              //8

    hila::out0 << "  config.initialCondition is "
	       << config.initialCondition
	       << "\n";
    
    config.variance_sigma = parameters.get("sigma");
    
    config.seed = parameters.get("seed");
    config.IniMod = parameters.get("IniMod");
    config.Inilc = parameters.get("Inilc");

    /*----------------------------------------*/
    /*       parameters update strategies     */
    /*----------------------------------------*/
    config.item = parameters.get_item("category",{"fixed", "computed", "interpolated"});
    if (config.item == 0){
      config.alpha = parameters.get("alpha");
      config.beta1 = parameters.get("beta1");
      config.beta2 = parameters.get("beta2");
      config.beta3 = parameters.get("beta3");
      config.beta4 = parameters.get("beta4");
      config.beta5 = parameters.get("beta5");
    }
    else if (config.item == 1) {
      //config.T = parameters.get("T");
      config.dT_from_TAB = parameters.get("dT_from_TAB");
      config.p = parameters.get("p");
      config.T = (MP.tAB_RWS(config.p) * MP.Tcp_mK(config.p)) + config.dT_from_TAB;
      hila::out0 << "Tab: "<< MP.tAB_RWS(config.p)
	         << ", p:" << config.p
	         << ", T:" << config.T
	         << ", dT_from_TAB:" << config.dT_from_TAB << "\n";
    }
    else {
      const std::string in_file = parameters.get("params_file"); 

      std::fstream params_stream;
      params_stream.open(in_file, std::ios::in);
      std::string line;

      while (std::getline(params_stream, line))
	{
	  std::stringstream ss(line);
	  real_t a, b, c;
	  if (ss >> a >> b >> c)
	    {
	      t_v.push_back(a);
	      T_v.push_back(b);
	      p_v.push_back(c);
	    }
	}
    }
    /*----------------------------------------*/
    /* parameters update strategies end here  */
    /*----------------------------------------*/   
    config.tStats = parameters.get("tStats");
    config.nOutputs = parameters.get("nOutputs");

    // output_file is the saving path of output file, which offered in congigration file
    const std::string output_file = parameters.get("output_file");

    // xdmf file name, which is provided through config file
    config.xmf2_fname = parameters.get("xmf2_file");

    /*if(config.positions==1)
      {
	config.npositionout = parameters.get("npositionout");
	config.write_phases = parameters.get_item("write_phases",{"no","yes"});
	config.write_eigen = parameters.get_item("write_eigen",{"no","yes"});
	}*/

    /*----------------------------------------*/
    /* >>>>>>>  boundary conditions  <<<<<<<<<*/
    /*----------------------------------------*/
    config.BCs1 = parameters.get_item("BCs1",{"periodic",
					       "AB",
					       "PairBreaking",
                                               "PB_y",
                                               "PairB_yz",
                                               "BB",
                                               "phaseVortices"});
    
    config.BCs2 = parameters.get_item("BCs2",{"periodic",
					      "AB",
					      "PairBreaking",
                                              "PB_y",
                                              "PairB_yz",
                                              "BB",
                                              "phaseVortices"});
    
    // config.Wn = parameters.get("BoundaryPhaseWindingNO");
    
    config.BCchangec = parameters.get("BCchangec");
    
    /*----------------------------------------*/
    /*       B-B domain wall BC data          */
    /*----------------------------------------*/
    // std::vector<real_t> tmp3 = parameters.get("BDiagMatrix_Left");
    // std::vector<real_t> tmp4 = parameters.get("BDiagMatrix_Right");    
    // if (tmp3.size() == 3 && tmp4.size() == 3)
    //   {
    // 	config.BLeft_11 = tmp3[0]; config.BLeft_22 = tmp3[1]; config.BLeft_33 = tmp3[2];
    // 	config.BRight_11 = tmp4[0]; config.BRight_22 = tmp4[1]; config.BRight_33 = tmp4[2];	
    //   }


    
    config.useTbath = parameters.get_item("useTbath",{"no","yes"});

    /*----------------------------------------*/
    /*       insitu randering parameters      */
    /*----------------------------------------*/
    // config.A_matrix_output             = parameters.get_item("A_matrix_output",{"no","yes"});
    // config.hdf5_A_matrix_output        = parameters.get_item("hdf5_A_matrix_output",{"no","yes"});
    // config.hdf5_trA_output             = parameters.get_item("hdf5_trA_output",{"no","yes"});
    // config.hdf5_eigvA_output           = parameters.get_item("hdf5_eigvA_output",{"no","yes"});
    // config.hdf5_mass_current_output    = parameters.get_item("hdf5_mass_current_output",{"no","yes"});
    // config.hdf5_spin_current_output    = parameters.get_item("hdf5_spin_current_output",{"no","yes"});        

    // config.do_gapA_clip         = parameters.get_item("do_gapA_clip",{"no","yes"});
    // config.do_gapA_isosurface   = parameters.get_item("do_gapA_isosurface",{"no","yes"});
    // config.do_gapA_3slice       = parameters.get_item("do_gapA_3slice",{"no","yes"});
    // config.do_fe_slice          = parameters.get_item("do_fe_slice",{"no","yes"});
    // config.do_gapA_slice        = parameters.get_item("do_gapA_slice",{"no","yes"});            
    
    config.clamp_bias_gapMin = parameters.get("clamp_bias_gapMin");
    config.clamp_bias_gapMax = parameters.get("clamp_bias_gapMax");
    // config.clamp_fed_Min = parameters.get("clamp_fed_Min");
    // config.clamp_fed_Max = parameters.get("clamp_fed_Max");    
    config.camera_azi = parameters.get("camera_azi");
    config.camera_ele = parameters.get("camera_ele");
    /*----------------------------------------*/
    /*  insitu randering parameters end       */
    /*----------------------------------------*/
    
    config.dt = config.dx * config.dtdxRatio;
    t = config.tStart;

    // setup the hila lattice geometry 
    CoordinateVector box_dimensions = {config.lx, config.ly, config.lz};
    lattice.setup(box_dimensions);
    hila::seed_random(config.seed);

    return output_file;
    
} // allocate() function ends here


void glsol::initialize() {

  Matep MP;
  real_t Tp[2];
  update_Tp(t, Tp);
  
  int Nx = config.lx;
  int Ny = config.ly;
  int Nz = config.lz;
  
  real_t dx = config.dx;

  real_t ttc = MP.Tcp_mK(Tp[1]);
  hila::out0 <<"T_AB: "<<MP.tAB_RWS(Tp[1])*ttc<<"\n";
  
  switch (config.initialCondition) {
    
  case 0: {
    pi = 0;                            
    real_t gap = MP.gap_B_td(Tp[1], Tp[0]);
    onsites(ALL) {                     
      A[X] = hila::gaussrand();
      A[X] = gap * A[X]/A[X].norm();   
    }

    hila::out0 << "Components randomly created \n";

    break;
    }
  case 1: {
    auto kA = A;
    real_t gap = config.IniMod;        //MP.gap_B_td(Tp[1], Tp[0]);
    real_t lc = config.Inilc;          //1.0/sqrt(abs(config.alpha));

    hila::out0 << "Correlation length in ICs: "<< lc <<"\n";
    
    onsites (ALL) {
            real_t constant = pow(gap, 2.0) * pow(2.0 * M_PI, 1.5) * pow(lc, 3.0)/(9.0 * Nx * Ny * Nz * dx * dx * dx);
            real_t kSqu;
            real_t std;
            kSqu = 0.0;
            auto k = X.coordinates();

            kSqu = pow(sin(M_PI * k.e(0) / Nx), 2.0)+pow(sin(M_PI * k.e(1) / Ny), 2.0)+pow(sin(M_PI * k.e(2) / Nz), 2.0); 
            kSqu *= pow(2.0 / dx, 2.0);

            if (kSqu > 0.0) {
                std = sqrt(0.5 * constant *
                           exp(-0.5 * kSqu * lc * lc));

		kA[X].gaussian_random(std);
		
            } else {
	      kA[X]=0;
	    }
    }

        FFT_field(kA, A, fft_direction::back);
	
	onsites (ALL)
	  {
	    // if (A[X].norm()>0.0){A[X]=gap*A[X]/A[X].norm();}
	    A[X]=A[X]/sqrt((A[X]*A[X].dagger()).trace());
	  }
	
        pi[ALL] = 0;

        hila::out0 << "k space generation \n";

        break;
  }

  case 2: {
    pi = 0.;                            
    onsites(ALL) {                     
      A[X] = sqrt(0.1) * hila::gaussrand();
    }

    hila::out0 << " normal-phase-real-1 created \n";
    break;    

  }
  case 3: {
    pi = 0.;
    onsites(ALL) {
      /*foralldir(al) foralldir(i){
        A[X].e(al,i).real().gaussian_random();
	}*/
      A[X].gaussian_random();
      A[X] = A[X].real();
    }
    hila::out0 << " normal-phase-real-2 created \n";
    break;    
  }

  case 4: {
    pi = 0.;
    onsites(ALL) {
      foralldir(al) foralldir(i){
	A[X].e(al,i) = sqrt(config.variance_sigma) * hila::gaussian_random<Complex<real_t>>();
      } // doralldir end here
    } // onsites(ALL) end here
  
    hila::out0 << " normal-phase-complex created \n";
    break;
      
  }
    
  case 5: {
    pi = 0;
    real_t gap = MP.gap_B_td(Tp[1], Tp[0]);
    hila::out0 <<"Gap B: "<<gap<<"\n";
    onsites(ALL) {
      foralldir(d1)foralldir(d2){

	if (d1==d2){
	  A[X].e(d1,d2).re = 1.0; //hila::gaussrand() hila::random()
	  A[X].e(d1,d2).im = 0.0;
	}
	else {
	  A[X].e(d1,d2).re = 0.0;
	  A[X].e(d1,d2).im = 0.0;}
      }
      A[X] = gap * A[X]/A[X].norm();
    }

    hila::out0 << "Pure B phase \n";

    break;
    }
    
    case 6: {
    pi = 0;
    real_t gap = MP.gap_A_td(Tp[1], Tp[0]);
    hila::out0<<"Gap A: "<<gap<<"\n";
    onsites(ALL) {

        A[X] = sqrt(config.variance_sigma) * hila::gaussrand();
	A[X].e(0,0).re = 1.0 + sqrt(config.variance_sigma) * hila::gaussrand();
	A[X].e(0,1).im = 1.0 + sqrt(config.variance_sigma) * hila::gaussrand();
      
        A[X] = gap * A[X]/sqrt(2.0);
    }

    hila::out0 << "Aphase_partial is created \n";

    break;
    }

  case 7: {
    pi = 0.;
    real_t gap = MP.gap_A_td(Tp[1], Tp[0]);
    hila::out0<<"Gap A: "<<gap<<"\n";
    
    onsites(ALL) {
      foralldir(al) foralldir(i){
	A[X].e(al,i) = sqrt(config.variance_sigma) * hila::gaussian_random<Complex<real_t>>();
	
	if ((al==0) && (i==0)) {
	  A[X].e(al,i).re=A[X].e(al,i).re + 1.;
	}
	else if ((al==0) && (i==1)) {
	  A[X].e(al,i).im=A[X].e(al,i).im + 1.;
	} // put bulk A-phase elements into random matrix

	A[X].e(al,i)=(A[X].e(al,i)/sqrt(2.)) * gap;	
      } // doralldir end here
    } // onsites(ALL) end here

    //A[ALL]=A[x].asArray()
    hila::out0 << "Aphase_full is created \n";

    break;
  } // case 7: Aphase_full

  case 8: {
    pi = 0;
    real_t gapa = MP.gap_A_td(Tp[1], Tp[0]);
    real_t gapb = MP.gap_B_td(Tp[1], Tp[0]);    
    // hila::out0<<"Gap A: "<<gap<<"\n";
    // if (X.coordinate(e_x) == 0 or X.coordinate(e_x) == 1)
    
    onsites (ALL) {
    if (
	(X.coordinate(e_x) <= 10 or X.coordinate(e_x) >= (config.lx - 10))
	&& ((X.coordinate(e_y) <= 74 and X.coordinate(e_y) >= 54))
	&& ((X.coordinate(e_z) <= 74 and X.coordinate(e_z) >= 54))
       )
      {	
	    foralldir(d1)foralldir(d2){
	      if (d1==d2){
		A[X].e(d1,d2).re = 1.0;
		A[X].e(d1,d2).im = 0.0;
	      }
	      else {
		A[X].e(d1,d2).re = 0.0;
		A[X].e(d1,d2).im = 0.0;}
              }
	    A[X] = gapb * A[X]/sqrt(3.0);
       }
    //else if (X.coordinate(e_z) == (config.lz - 1) or X.coordinate(e_z) == (config.lz - 2))
    else if (
	     (X.coordinate(e_x) > 10 and X.coordinate(e_x) < (config.lx - 10))
	     || (
    	         (X.coordinate(e_x) <= 10 or X.coordinate(e_x) >= (config.lx - 10))
		 && (!(
                        (X.coordinate(e_y) <= 74 and X.coordinate(e_y) >= 54)
  	                && (X.coordinate(e_z) <= 74 and X.coordinate(e_z) >= 54)
                      ))
                )
	    )    
	  {
	    foralldir(d1)foralldir(d2){
	      if (d1==0 && d2==0){
		A[X].e(d1,d2).re = 1.0;
		A[X].e(d1,d2).im = 0.0;
	      }
	      else if (d1==0 && d2==1){
		A[X].e(d1,d2).re = 0.0;  
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

    hila::out0 << "B domain in supercooling A \n";

    break;
    }
    
  default: {

    // #pragma hila ast_dump
    pi = 0.0; //set derivative matrix to zero
    onsites (ALL) {
      A[X].fill(1.0);
    }
    
    hila::out0 << "Field matrix set to 1 everywhere \n";

    break;
  }
  }

} // initialize() call end here

void glsol::update_Tp (real_t t, real_t Tp[2]) {

   
  if (config.item ==2)
    {
      int k1,k2;
      int size = t_v.size();

      if (t == t_v[size -1]) {
	k1 = size-1;
	k2 = size-1;}
      else {
	for (int i=1; i<size; i++) {
	  
	  if (t >= t_v[i-1] && t < t_v[i]) {
	    k1 = i-1;
	    k2 = i; }
	}
      }
      
      if (k1 == k2 && k1 == size-1){
	Tp[0] = T_v[size-1];
	Tp[1] = p_v[size-1];}
      else {
	Tp[0]= T_v[k1] + ((T_v[k2]-T_v[k1])/(t_v[k2]-t_v[k1])) * (t - t_v[k1]);
	Tp[1]= p_v[k1] + ((p_v[k2]-p_v[k1])/(t_v[k2]-t_v[k1])) * (t - t_v[k1]);}
    }
  else if (config.item == 1)
    {
      Tp[0] = config.T;
      Tp[1] = config.p;
    }
  else if (config.item == 0)
    {
      Tp[0] = 0.0;
      Tp[1] = 0.0;
    }

 
} // update_Tp function ends here

void glsol::update_params() {

  Matep MP;
  real_t Tp[2];
  real_t gapa,gapb;
  real_t fa,fb;
  
  if (config.item ==1 && t == config.tStart){
    tc = MP.Tcp_mK(config.p);
    gapa = MP.gap_A_td(config.p, config.T);
    gapb = MP.gap_B_td(config.p, config.T);
    fa = MP.f_A_td(config.p, config.T);
    fb = MP.f_B_td(config.p, config.T);
    
    hila::out0<<"T: "<<config.T<<" p:"<<config.p<<"\n";
    hila::out0<<"Tc: "<<tc<<"\n";
    hila::out0<<"Gap_A: "<<gapa<<" Gap_B: "<<gapb<< "\n";
    hila::out0<<"Bulk_A: "<<fa<<" Bulk_B: "<<fb<<"\n";
    
    config.alpha = MP.alpha_td(config.p, config.T);
    config.beta1 = MP.beta1_td(config.p, config.T);
    config.beta2 = MP.beta2_td(config.p, config.T);
    config.beta3 = MP.beta3_td(config.p, config.T);
    config.beta4 = MP.beta4_td(config.p, config.T);
    config.beta5 = MP.beta5_td(config.p, config.T);

    hila::out0 << config.alpha << " " << config.beta1 << " " << config.beta2 << " " << config.beta3 << " " << config.beta4 << " " << config.beta5<< "\n";

  }
  else if (config.item ==2){
    update_Tp(t, Tp);
    tc = MP.Tcp_mK(Tp[1]);
    gapa = MP.gap_A_td(Tp[1], Tp[0]);
    gapb = MP.gap_B_td(Tp[1], Tp[0]);
    fa = MP.f_A_td(Tp[1], Tp[0]);
    fb = MP.f_B_td(Tp[1], Tp[0]);

    hila::out0<<"t="<<t<<"\n";
    hila::out0<<"T: "<<Tp[0]<<" p:"<<Tp[1]<<"\n";
    hila::out0<<"Tc: "<<tc<<"\n";
    hila::out0<<"Gap_A: "<<gapa<<" Gap_B: "<<gapb<< "\n";
    hila::out0<<"Bulk_A: "<<fa<<" Bulk_B: "<<fb<<"\n";
    
    config.alpha = MP.alpha_td(Tp[1], Tp[0]);
    config.beta1 = MP.beta1_td(Tp[1], Tp[0]);
    config.beta2 = MP.beta2_td(Tp[1], Tp[0]);
    config.beta3 = MP.beta3_td(Tp[1], Tp[0]);
    config.beta4 = MP.beta4_td(Tp[1], Tp[0]);
    config.beta5 = MP.beta5_td(Tp[1], Tp[0]);
  }
} // update_params function ends here

void glsol::write_moduli() {

   // real_t a = scaleFactor(t);

  double Amod = 0.0/*, AGap = 0.0*/;
    double pimod = 0.0;
    
    real_t Tp[2];

    hila::set_allreduce(false);
    onsites (ALL) {
        Amod += A[X].norm();
        pimod += pi[X].norm();

    }
        
    update_Tp(t, Tp);
    
    if (hila::myrank() == 0) {
        config.stream << t << ", "
	              << tc << ", "
	              << Tp[0] << ", " << Tp[1] << ", "
	              << config.alpha << " " << config.beta1 << " " <<	config.beta2 << " " <<	config.beta3 << " " <<	config.beta4 << " "    <<	config.beta5 << " "
	    << Amod / lattice.volume() << " " << pimod / lattice.volume()
	              << std::endl;
    }
} // write_moduli ends at here

void glsol::write_energies() {

  Matep MP;
  
  Complex<double> suma(0),sumb2(0),sumb3(0),sumb4(0),sumb5(0);
  Complex<double> suma_we(0),sumb2_we(0),sumb3_we(0),sumb4_we(0),sumb5_we(0);
  Complex<double> sumk1(0),sumk2(0),sumk3(0);
  Complex<double> sumk1_we(0),sumk2_we(0),sumk3_we(0);
  real_t sumb1 = 0;
  Complex<double> sumb1_we = 0;
  Complex<double> sumkin(0);
  Complex<double> sumkin_we(0);
  real_t Tp[2];

  update_Tp(t, Tp);

  real_t ebfe=fmin(MP.f_A_td(Tp[1], Tp[0]),MP.f_B_td(Tp[1], Tp[0]));
  
    hila::set_allreduce(false);
    onsites(ALL) {

      Complex<double> a(0),b2(0),b3(0),b4(0),b5(0);
      Complex<double> kin(0);
      Complex<double> k1(0), k2(0), k3(0);
      Complex<double> bfe(0);
      double b1 = 0;
      
      a = config.alpha * (A[X]*A[X].dagger()).trace();
      b1 = config.beta1 * ((A[X]*A[X].transpose()).trace()).squarenorm();
      b2 = config.beta2 * ((A[X]*A[X].dagger()).trace()*(A[X]*A[X].dagger()).trace());
      b3 = config.beta3 * ((A[X]*A[X].transpose()*A[X].conj()*A[X].dagger()).trace());
      b4 = config.beta4 * ((A[X]*A[X].dagger()*A[X]*A[X].dagger()).trace());
      b5 = config.beta5 * ((A[X]*A[X].dagger()*A[X].conj()*A[X].transpose()).trace());

      bfe = a + b1 + b2 + b3 + b4 + b5 - ebfe;

      kin = (pi[X]*pi[X].dagger()).trace();
      
      foralldir(j) foralldir (k) foralldir(al){
	k1 += (A[X + k].column(j) - A[X - k].column(j)).e(al)
	      * (A[X + k].conj().column(j) - A[X - k].conj().column(j)).e(al)/(4.0*config.dx*config.dx);
	k2 += (A[X + j].column(j) - A[X - j].column(j)).e(al)
	      * (A[X + k].conj().column(k) - A[X - k].conj().column(k)).e(al)/(4.0*config.dx*config.dx);
	k3 += (A[X + k].column(j) - A[X - k].column(j)).e(al)
	      * (A[X + j].conj().column(k) - A[X - j].conj().column(k)).e(al)/(4.0*config.dx*config.dx);
      }

      sumkin += kin;
      sumkin_we += bfe * kin;

      sumk1 += k1;
      sumk1_we += bfe * k1;

      sumk2 += k2;
      sumk2_we += bfe *	k2;

      sumk3 += k3;
      sumk3_we += bfe *	k3;
      
      suma += a;
      suma_we += bfe * a;

      sumb1 += b1;
      sumb1_we += bfe * b1;

      sumb2 += b2;
      sumb2_we += bfe * b2;

      sumb3 += b3;
      sumb3_we += bfe * b3;

      sumb4 += b4;
      sumb4_we += bfe * b4;

      sumb5 += b5;
      sumb5_we += bfe * b5;
    }

      if (hila::myrank() == 0) {
        double vol = lattice.volume();
        config.stream << t << " "
	              << sumkin.re / vol << " " << sumkin.im / vol << " "
	              << sumkin_we.re / vol << " " << sumkin_we.im / vol << " "
		      << sumk1.re / vol << " " << sumk1.im / vol << " "
	              << sumk1_we.re / vol << " " << sumk1_we.im / vol << " "
		      << sumk2.re / vol	<< " " << sumk2.im / vol << " "
	              << sumk2_we.re / vol << " " << sumk2_we.im / vol << " "
		      << sumk3.re / vol	<< " " << sumk3.im / vol << " "
	              << sumk3_we.re / vol << " " << sumk3_we.im / vol << " "
		      << suma.re / vol << " " << suma.im / vol << " "
	              << suma_we.re / vol << " " << suma_we.im / vol << " "
		      << sumb1 / vol << " "
	              << sumb1_we.re / vol << " " << sumb1_we.im / vol << " "
		      << sumb2.re / vol << " " << sumb2.im / vol << " "
	              << sumb2_we.re / vol << " " << sumb2_we.im / vol << " "
		      << sumb3.re / vol << " " << sumb3.im / vol << " "
	              << sumb3_we.re / vol << " " << sumb3_we.im / vol << " "
		      << sumb4.re / vol << " " << sumb4.im / vol << " "
	              << sumb4_we.re / vol << " " << sumb4_we.im / vol << " "
		      << sumb5.re / vol << " " << sumb5.im / vol << " "
	              << sumb5_we.re / vol << " " << sumb5_we.im / vol << " "
	              << std::endl;
      }

       hila::out0 << "Energy done \n";

} // write_energies function ends here

void glsol::write_phases() {

  real_t ph0(0),ph1(0),ph2(0),ph3(0),ph4(0),ph5(0),ph6(0),ph7(0),ph8(0),ph9(0);

  hila::set_allreduce(false);

  onsites (ALL) {

    real_t R1,R2,R3,R4,R5;
    real_t p1,p2,p3,p4,p5,p6,p7,p8;
    real_t error=1.0/5.0;
    int phase;
    real_t p;
    phi_t Ac;

    if (real((A[X]*A[X].dagger()).trace()) > 0.0)
      {
	Ac=A[X]/sqrt((A[X]*A[X].dagger()).trace());;
      }
    else
      {
	Ac=A[X];
      }

    R1 = ((Ac*Ac.transpose()).trace()).squarenorm();

    R2 = real(((Ac*Ac.dagger()).trace()*(Ac*Ac.dagger()).trace()));

    R3 = real(((Ac*Ac.transpose()*Ac.conj()*Ac.dagger()).trace()));

    R4 = real(((Ac*Ac.dagger()*Ac*Ac.dagger()).trace()));

    R5 = real(((Ac*Ac.dagger()*Ac.conj()*Ac.transpose()).trace()));

    p1=abs(R1-1.0) + abs(R3-1.0/3.0) + abs(R4-1.0/3.0) + abs(R5-1.0/3.0);

    p2=abs(R1-1.0) + abs(R3-1.0/2.0) + abs(R4-1.0/2.0) + abs(R5-1.0/2.0);

    p3=abs(R1-1.0) + abs(R3-1.0) + abs(R4-1.0) + abs(R5-1.0);

    p4=abs(R1-0.0) + abs(R3-1.0/3.0) + abs(R4-1.0/3.0) + abs(R5-1.0/3.0);

    p5=abs(R1-0.0) + abs(R3-1.0/2.0) + abs(R4-1.0/2.0) + abs(R5-1.0/2.0);

    p6=abs(R1-0.0) + abs(R3-0.0) + abs(R4-1.0) + abs(R5-1.0);

    p7=abs(R1-0.0) + abs(R3-1.0) + abs(R4-1.0) + abs(R5-0.0);

    p8=abs(R1-0.0) + abs(R3-0.0) + abs(R4-1.0) + abs(R5-0.0);

    if(p1<p2 and p1<p3 and p1<p4 and p1<p5 and p1<p6 and p1<p7 and p1<p8 and p1<error)
      {
	ph1+=1;
      }
    else if (p2<p1 and p2<p3 and p2<p4 and p2<p5 and p2<p6 and p2<p7 and p2<p8 and p2<error)
      {
	ph2+=1;
      }
    else if (p3<p2 and p3<p1 and p3<p4 and p3<p5 and p3<p6 and p3<p7 and p3<p8 and p3<error)
      {
	ph3+=1;
      }
    else if (p4<p2 and p4<p3 and p4<p1 and p4<p5 and p4<p6 and p4<p7 and p4<p8 and p4<error)
      {
	ph4+=1;
      }
    else if (p5<p2 and p5<p3 and p5<p4 and p5<p1 and p5<p6 and p5<p7 and p5<p8 and p5<error)
      {
        ph5+=1;
      }
    else if (p6<p2 and p6<p3 and p6<p4 and p6<p5 and p6<p1 and p6<p7 and p6<p8 and p6<error)
      {
	ph6+=1;
      }
    else if (p7<p2 and p7<p3 and p7<p4 and p7<p5 and p7<p6 and p7<p1 and p7<p8 and p7<error)
      {
	ph7+=1;
      }
    else if (p8<p2 and p8<p3 and p8<p4 and p8<p5 and p8<p6 and p8<p7 and p8<p1 and p8<error)
      {
	ph8+=1;
      }
    else
      {
	ph0+=1;
      }
  }

  if (hila::myrank() == 0)
    {
      double vol = lattice.volume();
      config.stream << ph0/vol << " " << ph1/vol << " " << ph2/vol << " " << ph3/vol << " " << ph4/vol << " " << ph5/vol << " " << ph6/vol << " " << ph7/vol << " " << ph8/vol << " " << ph9/vol << "\n";
    }

} // write_phases() function ends here

void glsol::write_positions() {

  Matep MP;
  real_t Tp[2];
  update_Tp(t, Tp);

  real_t gapa = MP.gap_A_td(Tp[1], Tp[0]);
  real_t gapb = MP.gap_B_td(Tp[1], Tp[0]);
  
  hila::set_allreduce(false);

  if (config.write_phases==1)
    {
      Field<Vector<8,float>> data;
      
      onsites (ALL) {

	real_t R1,R2,R3,R4,R5;
	phi_t Ac;

	if (real((A[X]*A[X].dagger()).trace()) > 0.0)
	  {
	    Ac=A[X]/sqrt((A[X]*A[X].dagger()).trace());;
	  }
	else
	  {
	    Ac=A[X];
	  }
	    
	R1 = ((Ac*Ac.transpose()).trace()).squarenorm();

	R2 = real(((Ac*Ac.dagger()).trace()*(Ac*Ac.dagger()).trace()));

	R3 = real(((Ac*Ac.transpose()*Ac.conj()*Ac.dagger()).trace()));

	R4 = real(((Ac*Ac.dagger()*Ac*Ac.dagger()).trace()));

	R5 = real(((Ac*Ac.dagger()*Ac.conj()*Ac.transpose()).trace()));

	data[X].e(0)=abs(R1-1.0) + abs(R3-1.0/3.0) + abs(R4-1.0/3.0) + abs(R5-1.0/3.0);

	data[X].e(1)=abs(R1-1.0) + abs(R3-1.0/2.0) + abs(R4-1.0/2.0) + abs(R5-1.0/2.0);

	data[X].e(2)=abs(R1-1.0) + abs(R3-1.0) + abs(R4-1.0) + abs(R5-1.0);

	data[X].e(3)=abs(R1-0.0) + abs(R3-1.0/3.0) + abs(R4-1.0/3.0) + abs(R5-1.0/3.0);

	data[X].e(4)=abs(R1-0.0) + abs(R3-1.0/2.0) + abs(R4-1.0/2.0) + abs(R5-1.0/2.0);

	data[X].e(5)=abs(R1-0.0) + abs(R3-0.0) + abs(R4-1.0) + abs(R5-1.0);

	data[X].e(6)=abs(R1-0.0) + abs(R3-1.0) + abs(R4-1.0) + abs(R5-0.0);

	data[X].e(7)=abs(R1-0.0) + abs(R3-0.0) + abs(R4-1.0) + abs(R5-0.0);

      }

      data.write("points/phase-distances-t"+std::to_string(int(t/config.dt)),false);

    }

  if(config.write_eigen==1)
    {

      Field<Vector<3,double>> eval;
      Field<Matrix<3,3,Complex<double>>> evec;

      onsites(ALL){

	(A[X].dagger()*A[X]).eigen_jacobi(eval[X],evec[X],hila::sort::ascending);

      }

      eval.write("points/eigenvalues-t"+std::to_string(int(t/config.dt)),false);
      evec.write("points/eigenvectors-t"+std::to_string(int(t/config.dt)),false);
      
    }
      
} // write_positions() function ends here

void glsol::write_xdmf(){

    unsigned int rank_no = 0/*hila::myrank()*/;
    std::fstream xml_file;
    config.xdmf_out.open(config.xmf2_fname, std::ios::out);

    config.xdmf_out << "<?xml version=\"1.0\" ?>\n"
                    << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
                    << "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.0\">\n"
                    << "\n"
                    << "<Domain>\n"
                    << "<Grid GridType=\"Collection\" CollectionType=\"Collection\">"
                    << "\n";
    
    while(rank_no < hila::number_of_nodes()){
      const std::string fname = "rank_xmls/" + config.xmf2_fname + "_" + std::to_string(rank_no) + ".xml";  
      xml_file.open(fname, std::ios::in);
      config.xdmf_out << xml_file.rdbuf() << "\n" << std::endl;
      xml_file.close();
      
      ++rank_no;
    }

    config.xdmf_out << "</Grid>"
                    << "\n"
                    << "</Domain>\n"
                    << "</Xdmf>"
                    << std::endl;
    config.xdmf_out.close();
} // write_xdmf() function ends here

// void glsol::write_A_matrix_positions() {

//   std::ofstream stream_out;
  
//   const std::string fname = "A_matrix_output/t"+std::to_string(int(t/config.dt))+".dat";
//   stream_out.open(fname, std::ios::out);

//   hila::set_allreduce(false);
//   A.write(stream_out,false,8);
  
//   stream_out.close();

// }

void glsol::next() {

  Matep MP;
  real_t Tp[2];
  update_Tp(t, Tp);
  
  static hila::timer next_timer("timestep");
  Field<phi_t> deltaPi;
  Field<Vector<3,Complex<real_t>>> djAaj;

  real_t gapa = MP.gap_A_td(Tp[1], Tp[0]);
  real_t gapb = MP.gap_B_td(Tp[1], Tp[0]);

  int bc=config.boundaryConditions;
  //std::string bc=config.boundaryConditions;
  
  next_timer.start();

  /* --------------------------------------------------------------- */  
  /* >>>>>>>>>>  boundary condition handling block starts <<<<<<<<<< */
  /* --------------------------------------------------------------- */
  onsites (ALL) {

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
		A[X].e(d1,d2).im = 0.0;}
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
	       if (X.coordinate(e_x) == 0 or X.coordinate(e_x) == (config.lx - 1) or
	           X.coordinate(e_x) == 1 or X.coordinate(e_x) == (config.lx - 2) or
	           X.coordinate(e_y) == 0 or X.coordinate(e_y) == (config.ly - 1) or
	           X.coordinate(e_y) == 1 or X.coordinate(e_y) == (config.ly - 2) or
	           X.coordinate(e_z) == 0 or X.coordinate(e_z) == (config.lz - 1) or
	           X.coordinate(e_z) == 1 or X.coordinate(e_z) == (config.lz - 2))
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
      
    auto AxAt = A[X]*A[X].transpose();
    auto AxAd = A[X]*A[X].dagger();
    
    deltaPi[X] = - config.alpha*A[X]
      - 2.0*config.beta1*A[X].conj()*AxAt.trace() 
      - 2.0*config.beta2*A[X]*AxAd.trace()
      - 2.0*config.beta3*AxAt*A[X].conj()
      - 2.0*config.beta4*AxAd*A[X] 
      - 2.0*config.beta5*A[X].conj()*A[X].transpose()*A[X];

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

  onsites (ALL) {
    deltaPi[X] +=  (1.0 / (config.dx * config.dx)) * (A[X + e_x] + A[X - e_x]
						      + A[X + e_y] + A[X - e_y]
						      + A[X + e_z] + A[X - e_z]
						      - 6.0 * A[X]);
  } // DjDjAali term summation

    //onsites (ALL) {deltaPi[X] *= config.dt;} // I think that this is the problem, multiplication with respect to dt

  if (t < config.tdif)
    {
      pi[ALL] = deltaPi[X]/(config.difFac);
      t += config.dt/config.difFac;
    }
  else if (t < config.tdis && /*config.gamma*/ config.gamma.real() >= 0 )
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

void glsol::next_bath() {

  Matep MP;
  real_t Tp[2];
  update_Tp(t, Tp);

  static hila::timer next_timer("timestep");
  Field<phi_t> deltaPi;
  Field<Vector<3,Complex<real_t>>> djAaj;

  real_t gapa = MP.gap_A_td(Tp[1], Tp[0]);
  real_t gapb = MP.gap_B_td(Tp[1], Tp[0]);

  real_t tb = Tp[0]/tc;
  //real_t sig = sqrt(2.0*tb*config.gamma);
  Complex<real_t> sig = sqrt(2.0*tb*config.gamma);  
  //real_t ep2 = 1.0-exp(-2.0*config.gamma*config.dt) ;
  Complex<real_t> ep2 = 1.0-exp(-2.0*config.gamma*config.dt) ;

  int bc=config.boundaryConditions;
  //std::string bc=config.boundaryConditions;

  next_timer.start();

  onsites (ALL) {

    A[X] += config.dt * pi[X];

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
                A[X].e(d1,d2).im = 0.0;}
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
                A[X].e(d1,d2).re = 0.0;
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
        if (X.coordinate(e_x) == 0 or X.coordinate(e_x) == (config.lx - 1) or
            X.coordinate(e_x) == 1 or X.coordinate(e_x) == (config.lx - 2) or
            X.coordinate(e_y) == 0 or X.coordinate(e_y) == (config.ly - 1) or
            X.coordinate(e_y) == 1 or X.coordinate(e_y) == (config.ly - 2) or
            X.coordinate(e_z) == 0 or X.coordinate(e_z) == (config.lz - 1) or
            X.coordinate(e_z) == 1 or X.coordinate(e_z) == (config.lz - 2))
          {
            A[X]=0.0;
          }
      }
  }

  onsites (ALL) {

    auto AxAt = A[X]*A[X].transpose();
    auto AxAd = A[X]*A[X].dagger();

    deltaPi[X] = - config.alpha*A[X]
      - 2.0*config.beta1*A[X].conj()*AxAt.trace()
      - 2.0*config.beta2*A[X]*AxAd.trace()
      - 2.0*config.beta3*AxAt*A[X].conj()
      - 2.0*config.beta4*AxAd*A[X]
      - 2.0*config.beta5*A[X].conj()*A[X].transpose()*A[X];

  }

  onsites(ALL) {
    djAaj[X] = 0;
    foralldir(j) {
      djAaj[X] += A[X + j].column(j) - A[X - j].column(j);
    }
  }

  onsites(ALL) {
    phi_t mat;
    foralldir(d) {
      auto col = djAaj[X+d] - djAaj[X-d];
      for (int i=0; i<NDIM; i++) mat.e(i,d) = col[i];
    }

    deltaPi[X] += (1.0/(2.0*(config.dx*config.dx)))*mat;
  }

  onsites (ALL) {
    deltaPi[X] += (1.0/(4.0*config.dx*config.dx)) * (A[X + e_x] + A[X - e_x]
                                                     + A[X + e_y] + A[X - e_y]
                                                     + A[X + e_z] + A[X - e_z]
                                                     - 6.0*A[X]);
     deltaPi[X] += hila::gaussrand()*sig;

  }

    //onsites (ALL) {deltaPi[X] *= config.dt;} // I think that this is the problem, multiplication with respect to dt                                                                                                                                                                                         
  if (t < config.tdif)
    {
      pi[ALL] = deltaPi[X]/(config.difFac);
      t += config.dt/config.difFac;
    }
  else if (t < config.tdis && /*config.gamma*/config.gamma.real() > 0 )
    {
      pi[ALL] = pi[X] + (deltaPi[X] - 2.0 * config.gamma * pi[X])*(config.dt/2.0);

      pi[ALL] = sqrt(1.0-ep2)*pi[X] + sqrt(ep2)*hila::gaussrand();


      deltaPi[ALL] = 0.0;

      onsites (ALL) {

        auto AxAt = A[X]*A[X].transpose();
        auto AxAd = A[X]*A[X].dagger();

        deltaPi[X] = - config.alpha*A[X]
          - 2.0*config.beta1*A[X].conj()*AxAt.trace()
          - 2.0*config.beta2*A[X]*AxAd.trace()
          - 2.0*config.beta3*AxAt*A[X].conj()
          - 2.0*config.beta4*AxAd*A[X]
          - 2.0*config.beta5*A[X].conj()*A[X].transpose()*A[X];
      }

      onsites(ALL) {
        djAaj[X] = 0;
        foralldir(j) {
          djAaj[X] += A[X + j].column(j) - A[X - j].column(j);
        }
      }

      onsites(ALL) {
        phi_t mat;
        foralldir(d) {
          auto col = djAaj[X+d] - djAaj[X-d];
          for (int i=0; i<NDIM; i++) mat.e(i,d) = col[i];
        }

        deltaPi[X] += (1.0/(2.0*(config.dx*config.dx)))*mat;
      }

      onsites (ALL) {
        deltaPi[X] += (1.0/(4.0*config.dx*config.dx)) * (A[X + e_x] + A[X - e_x]
                                                     + A[X + e_y] + A[X - e_y]
                                                     + A[X + e_z] + A[X - e_z]
							 - 6.0*A[X]);
        deltaPi[X] += hila::gaussrand()*sig;
      }


      pi[ALL] = pi[X] + (deltaPi[X] - 2.0 * config.gamma * pi[X])*(config.dt/2.0);

      t += config.dt;
    }
  else
    {
      pi[ALL] = pi[X] + deltaPi[X]*config.dt;
      t += config.dt;
    }


  next_timer.stop();

} // next_bath() ends here

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
#include "matep.hpp"

/*----------------------------------------------------------------------*/
/***  include Ascent & Canduit macro for in situ rank rendering      ****/
#if defined USE_ASCENT

#include "ascent.hpp"
#include "conduit_blueprint.hpp"

#endif
/*----------------------------------------------------------------------*/

// Definition of the field that we will use
using real_t = float;                          // or double ?
using phi_t = Matrix<3,3,Complex<real_t>>;     // saves the trouble of writing this every time

// Container for simulation parameters and methods
class he3sim{

public:
  he3sim() = default;                   // default constructor
  //he3sim(): matep {};                   // constructor for PLTS-2000 scale
  
  // read configration file and initiate scaling_sim.config 
  const std::string allocate(const std::string &fname, int argc, char **argv);
  
  void initialize();
  void update_params();
  void update_Tp (real_t t, real_t Tp[2]);

  void next();
  void next_bath();
  
  void write_moduli();
  void write_energies();
  void write_positions();
  void write_phases();
  void write_xdmf();
  
  void write_A_matrix_positions();             // output A-matrix after certain time interval
  void write_gapA_rankOut();
  void latticeCoordinate_output();             // lattice coordinates output with same sequence of A.write()

#if defined USE_ASCENT
    // In-Situ via Ascent
    void insitu_createMesh();
    void insitu_defineActions();
    void insitu_hdf5xdmf();
    void insitu_initialize();  
    void insitu_execute();
    void insitu_close();

    /*------------------------------*/
    Field<real_t> gapA;
    Field<real_t> feDensity;
    Field<real_t> u11, u12, u13, u21, u22, u23, u31, u32, u33;
    Field<real_t> v11, v12, v13, v21, v22, v23, v31, v32, v33;
   
    std::vector<real_t> gapAOrdered;
    std::vector<real_t> feDensityOrdered;
    std::vector<real_t> u11Ordered, u12Ordered, u13Ordered, u21Ordered, u22Ordered, u23Ordered, u31Ordered, u32Ordered, u33Ordered;
    std::vector<real_t> v11Ordered, v12Ordered, v13Ordered, v21Ordered, v22Ordered, v23Ordered, v31Ordered, v32Ordered, v33Ordered;
    /*------------------------------*/
  
    long long ghostVolume;
    long long latticeVolumeWithGhost;
    long long latticeVolume;
    unsigned char *ghostCellsMask;
    ascent::Ascent insitu;
    conduit::Node ascent_options;
    conduit::Node actions;
    conduit::Node mesh;
#endif
  
  Field<phi_t> A;
  Field<phi_t> pi;

  Matep matep;
  real_t t;
  real_t tc = 0;

  std::vector<real_t> t_v;
  std::vector<real_t> T_v;
  std::vector<real_t> p_v;
  
    struct config {
      int lx;
      int ly;
      int lz;
      real_t dx;
      real_t dt;
      real_t dtdxRatio;
      
      real_t tStart;
      real_t tEnd;

      real_t tdif;
      real_t difFac;
      real_t tdis;
      real_t gamma;

      int initialCondition;
      real_t variance_sigma;
      
      int seed;
      real_t IniMod;
      real_t Inilc;

      int item;
      real_t T;
      real_t dT_from_TAB;
      real_t p;
      real_t alpha;
      real_t beta1;
      real_t beta2;
      real_t beta3;
      real_t beta4;
      real_t beta5;

      real_t tStats;
      real_t nOutputs;

      std::fstream stream;
      std::fstream rank_Output;
      
      // xdmf file fstream
      std::fstream xdmf_out;
      std::string xmf2_fname;;
      int rank_stream_switch;

      int positions;
      int npositionout;
      int boundaryConditions;
      int useTbath;

      /* >>>>>>> insitu contral parameter <<<<<<<< */
      real_t clamp_bias_gapMin, clamp_bias_gapMax;
      real_t camera_azi, camera_ele;
      //unsigned int number_of_hdf;
      
      /* *******************************************/
      
      int write_phases;
      int write_eigen;
    } config;
};

const std::string he3sim::allocate(const std::string &fname, int argc, char **argv)
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
    config.gamma = parameters.get("gamma");
    config.initialCondition = parameters.get_item("initialCondition",{"gaussrand"
								      ,"kgaussrand"
								      ,"normal_phase"
								      ,"Bphase"
								      ,"Aphase"
                                                                      ,"BinA"});
    config.variance_sigma = parameters.get("sigma");
    
    config.seed = parameters.get("seed");
    config.IniMod = parameters.get("IniMod");
    config.Inilc = parameters.get("Inilc");

    // ******************************************************************************** //
    // >>>>>>>>>>>>>>>       different way for updating parameters      <<<<<<<<<<<<<<< //
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
      const std::string in_file = parameters.get("params_file"); // where is params_file ???

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
    // ******************************************************************************** //
    
    config.tStats = parameters.get("tStats");
    config.nOutputs = parameters.get("nOutputs");

    // output_file is the saving path of output file, which offered in congigration file
    const std::string output_file = parameters.get("output_file");

    // xdmf file name, which is provided through config file
    config.xmf2_fname = parameters.get("xmf2_file");

    config.rank_stream_switch = parameters.get_item("rank_output_stream",{"no", "yes"});
    /*if(config.positions==1)
      {
	config.npositionout = parameters.get("npositionout");
	config.write_phases = parameters.get_item("write_phases",{"no","yes"});
	config.write_eigen = parameters.get_item("write_eigen",{"no","yes"});
	}*/
    config.boundaryConditions = parameters.get_item("boundaryConditions",{"periodic",
									  "AB",
									  "PairBreaking",
                                                                          "PB_y",
                                                                          "PairB_yz"});
    config.useTbath = parameters.get_item("useTbath",{"no","yes"});

    /********************************************************************************/
    /* >>>>>>>>>>>>>          insitu randering parameters             <<<<<<<<<<<<<<*/
    config.clamp_bias_gapMin = parameters.get("clamp_bias_gapMin");
    config.clamp_bias_gapMax = parameters.get("clamp_bias_gapMax");
    config.camera_azi = parameters.get("camera_azi");
    config.camera_ele = parameters.get("camera_ele");
    //config.number_of_hdf = parameters.get("number_of_hdf");
    
    /**********************************************************************************/
    
    config.dt = config.dx * config.dtdxRatio;
    t = config.tStart;

    // setup the hila lattice geometry 
    CoordinateVector box_dimensions = {config.lx, config.ly, config.lz};
    lattice.setup(box_dimensions);
    hila::seed_random(config.seed);


    return output_file;
}


void he3sim::initialize() {

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
    pi = 0;                            
    //real_t gap = MP.gap_B_td(Tp[1], Tp[0]);
    onsites(ALL) {                     
      A[X] = sqrt(0.1) * hila::gaussrand();
    }

    hila::out0 << " nNormal phase created \n";

    break;    

  }    
  case 3: {
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
  case 4: {
    pi = 0;
    real_t gap = MP.gap_A_td(Tp[1], Tp[0]);
    hila::out0<<"Gap A: "<<gap<<"\n";
    onsites(ALL) {

        A[X] = sqrt(config.variance_sigma) * hila::gaussrand();
	A[X].e(0,0).re = 1.0 + sqrt(config.variance_sigma) * hila::gaussrand();
	A[X].e(0,1).im = 1.0 + sqrt(config.variance_sigma) * hila::gaussrand();
      
        A[X] = gap * A[X]/sqrt(2.0);
    }

    hila::out0 << "Pure A phase \n";

    break;
    }

  case 5: {
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
}

void he3sim::update_Tp (real_t t, real_t Tp[2]) {

   
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

 
}
void he3sim::update_params() {

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
}

void he3sim::write_moduli() {

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
}

void he3sim::write_energies() {

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
        config.stream << sumkin.re / vol << " " << sumkin.im / vol << " "
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
	              << sumb5_we.re / vol << " " << sumb5_we.im / vol << " ";
    }

       hila::out0 << "Energy done \n";

}

void he3sim::write_phases() {

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

}

void he3sim::write_positions() {

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
      
}

void he3sim::write_xdmf(){

    unsigned int rank_no = 0/*hila::myrank()*/;
    std::ifstream xml_file;
    config.xdmf_out.open(config.xmf2_fname, std::ios::out);

    config.xdmf_out << "<?xml version=\"1.0\" ?>\n"
                    << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
                    << "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.0\">\n"
                    << "\n"
                    << "<Domain>\n"
                    << "<Grid GridType=\"Collection\" CollectionType=\"Collection\">"
                    << "\n";
    while(rank_no < hila::number_of_nodes()){

      xml_file.open("rank_xmls/" + config.xmf2_fname + "_" + std::to_string(rank_no) + ".xml");
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
}

void he3sim::write_A_matrix_positions() {

  /*Matep MP;
    real_t Tp[2];
    update_Tp(t, Tp);*/

  std::ofstream stream_out;
  
  const std::string fname = "A_matrix_output/t"+std::to_string(int(t/config.dt))+".dat";
  stream_out.open(fname, std::ios::out);

  hila::set_allreduce(false);
  A.write(stream_out,false,8);
  
  stream_out.close();

  // ******************************************************************************** //

}

void he3sim::write_gapA_rankOut()
{
  onsites(ALL){
	if(config.rank_stream_switch == 1)
	  {
	    config.rank_Output << t << ", "
			       << X.coordinate(e_x) << ", " << X.coordinate(e_y) << ", " << X.coordinate(e_z) << ", "
	                       << sqrt((A[X]*A[X].dagger()).trace())
	                       << std::endl;
	  }
  }	
}  

void he3sim::latticeCoordinate_output()
{
  std::ofstream stream_out;
  
  const std::string fname_coordinates = "A_matrix_output/t"+std::to_string(int(t/config.dt))+"coordinates"+".dat";

  const long Nx = config.lx, Ny = config.ly, Nz = config.lz;
  const long NxNyNz = (Nx * Ny) * Nz;
  
  std::vector<long> lineNumList;

  for (long index = 0; index != NxNyNz; ++index)
    {
      lineNumList.push_back(index);
    }
  
  stream_out.open(fname_coordinates, std::ios::out);

  for (auto &i : lineNumList)
    {
      stream_out << " " << i % (Nx) << " " << (i/Nx) % Ny << " " << i/(Nx * Ny) << "\n";
    }

  stream_out.close();

}

/*--------------------------------------------------------------------------------*/
/*   >>>>>>>>>>>>>>>>>>    in-situ rank render functrions       <<<<<<<<<<<<<<<<< */
/*--------------------------------------------------------------------------------*/

#if defined USE_ASCENT
void he3sim::insitu_createMesh() {

    // Create a 3D mesh defined on a uniform grid of points
    // with a single vertex associated field named `gapA`
   
    // conduit::Node mesh;
    mesh["state/time"].set_external(&t);
    // mesh["state/cycle"].set_external(&step);
#if defined USE_MPI
    mesh["state/domain_id"] = lattice.mynode.rank;
#endif
    mesh["state/software"] = "He3Sim";
    mesh["state/title"] = "Bulk TDGL equation simulator";
    mesh["state/info"] = "In Situ rendering of order parameter field from He3Sim";

    // create the coordinate set
    mesh["coordsets/coords/type"] = "uniform";
    mesh["coordsets/coords/dims/i"] = lattice.mynode.size[0] + 2;
    mesh["coordsets/coords/dims/j"] = lattice.mynode.size[1] + 2;
    mesh["coordsets/coords/dims/k"] = lattice.mynode.size[2] + 2;

    // add origin and spacing to the coordset (optional)
    mesh["coordsets/coords/origin/x"] = ((lattice.mynode.min[0] - 1) * config.dx);
    mesh["coordsets/coords/origin/y"] = ((lattice.mynode.min[1] - 1) * config.dx);
    mesh["coordsets/coords/origin/z"] = ((lattice.mynode.min[2] - 1) * config.dx);
    mesh["coordsets/coords/spacing/dx"] = config.dx;
    mesh["coordsets/coords/spacing/dy"] = config.dx;
    mesh["coordsets/coords/spacing/dz"] = config.dx;
    
    // add the topology
    // this case is simple b/c it's implicitly derived from the coordinate set
    mesh["topologies/topo/type"] = "uniform";
    // reference the coordinate set by name
    mesh["topologies/topo/coordset"] = "coords";

    // create an vertex associated field named gapAOrdered
    mesh["fields/gapAOrdered/association"] = "vertex";
    mesh["fields/gapAOrdered/topology"] = "topo";
    mesh["fields/gapAOrdered/values"].set_external(gapAOrdered.data(), latticeVolumeWithGhost);

    // create an vertex associated field named feDensityOrdered
    mesh["fields/feDensityOrdered/association"] = "vertex";
    mesh["fields/feDensityOrdered/topology"] = "topo";
    mesh["fields/feDensityOrdered/values"].set_external(feDensityOrdered.data(), latticeVolumeWithGhost);

    /*----------------------------------------------------------------------*/
    /*---create an vertex associated field named of uxxOrdered vxxOrdered --*/
    mesh["fields/u11Ordered/association"] = "vertex";
    mesh["fields/u11Ordered/topology"] = "topo";
    mesh["fields/u11Ordered/values"].set_external(u11Ordered.data(), latticeVolumeWithGhost);

    mesh["fields/u12Ordered/association"] = "vertex";
    mesh["fields/u12Ordered/topology"] = "topo";
    mesh["fields/u12Ordered/values"].set_external(u12Ordered.data(), latticeVolumeWithGhost);

    mesh["fields/u13Ordered/association"] = "vertex";
    mesh["fields/u13Ordered/topology"] = "topo";
    mesh["fields/u13Ordered/values"].set_external(u13Ordered.data(), latticeVolumeWithGhost);

    mesh["fields/u21Ordered/association"] = "vertex";
    mesh["fields/u21Ordered/topology"] = "topo";
    mesh["fields/u21Ordered/values"].set_external(u21Ordered.data(), latticeVolumeWithGhost);

    mesh["fields/u22Ordered/association"] = "vertex";
    mesh["fields/u22Ordered/topology"] = "topo";
    mesh["fields/u22Ordered/values"].set_external(u22Ordered.data(), latticeVolumeWithGhost);

    mesh["fields/u23Ordered/association"] = "vertex";
    mesh["fields/u23Ordered/topology"] = "topo";
    mesh["fields/u23Ordered/values"].set_external(u23Ordered.data(), latticeVolumeWithGhost);

    mesh["fields/u31Ordered/association"] = "vertex";
    mesh["fields/u31Ordered/topology"] = "topo";
    mesh["fields/u31Ordered/values"].set_external(u31Ordered.data(), latticeVolumeWithGhost);

    mesh["fields/u32Ordered/association"] = "vertex";
    mesh["fields/u32Ordered/topology"] = "topo";
    mesh["fields/u32Ordered/values"].set_external(u32Ordered.data(), latticeVolumeWithGhost);

    mesh["fields/u33Ordered/association"] = "vertex";
    mesh["fields/u33Ordered/topology"] = "topo";
    mesh["fields/u33Ordered/values"].set_external(u33Ordered.data(), latticeVolumeWithGhost);

    mesh["fields/v11Ordered/association"] = "vertex";
    mesh["fields/v11Ordered/topology"] = "topo";
    mesh["fields/v11Ordered/values"].set_external(v11Ordered.data(), latticeVolumeWithGhost);

    mesh["fields/v12Ordered/association"] = "vertex";
    mesh["fields/v12Ordered/topology"] = "topo";
    mesh["fields/v12Ordered/values"].set_external(v12Ordered.data(), latticeVolumeWithGhost);

    mesh["fields/v13Ordered/association"] = "vertex";
    mesh["fields/v13Ordered/topology"] = "topo";
    mesh["fields/v13Ordered/values"].set_external(v13Ordered.data(), latticeVolumeWithGhost);

    mesh["fields/v21Ordered/association"] = "vertex";
    mesh["fields/v21Ordered/topology"] = "topo";
    mesh["fields/v21Ordered/values"].set_external(v21Ordered.data(), latticeVolumeWithGhost);

    mesh["fields/v22Ordered/association"] = "vertex";
    mesh["fields/v22Ordered/topology"] = "topo";
    mesh["fields/v22Ordered/values"].set_external(v22Ordered.data(), latticeVolumeWithGhost);

    mesh["fields/v23Ordered/association"] = "vertex";
    mesh["fields/v23Ordered/topology"] = "topo";
    mesh["fields/v23Ordered/values"].set_external(v23Ordered.data(), latticeVolumeWithGhost);

    mesh["fields/v31Ordered/association"] = "vertex";
    mesh["fields/v31Ordered/topology"] = "topo";
    mesh["fields/v31Ordered/values"].set_external(v31Ordered.data(), latticeVolumeWithGhost);

    mesh["fields/v32Ordered/association"] = "vertex";
    mesh["fields/v32Ordered/topology"] = "topo";
    mesh["fields/v32Ordered/values"].set_external(v32Ordered.data(), latticeVolumeWithGhost);

    mesh["fields/v33Ordered/association"] = "vertex";
    mesh["fields/v33Ordered/topology"] = "topo";
    mesh["fields/v33Ordered/values"].set_external(v33Ordered.data(), latticeVolumeWithGhost);        
    /*----------------------------------------------------------------------*/
    
    // // create an element associated field named ghostCells
    mesh["fields/ascent_ghosts/association"] = "element";
    mesh["fields/ascent_ghosts/topology"] = "topo";
    mesh["fields/ascent_ghosts/values"].set_external(ghostCellsMask, ghostVolume);

    // make sure the mesh we created conforms to the blueprint
    conduit::Node verify_info;
    if (!conduit::blueprint::mesh::verify(mesh, verify_info)) {
        hila::out0 << "Mesh Verify failed!\n";
        hila::out0 << verify_info.to_yaml() << '\n';
    } else {
        hila::out0 << "Mesh verify success!\n";
    }
};

void he3sim::insitu_defineActions() {

    // setup actions
    // conduit::Node actions;
    conduit::Node &add_act = actions.append();
    //add_act["action"] = "add_pipelines";
    //conduit::Node &pipelines = add_act["pipelines"];
    add_act["action"] = "add_scenes";
    conduit::Node &scenes = add_act["scenes"];

    /* >>>>>>>>>>>>>>   1st scene    <<<<<<<<<<<<< */
    
    // our first scene (named 's1') will render the field 'gapAOrdered'
    // to the file out_scene_ex1_render_var1.png
    scenes["s1/plots/p1/type"] = "pseudocolor";
    scenes["s1/plots/p1/field"] = "gapAOrdered";

    // color map clamping. min_value will be set to 0.0 if initialCondtion is 2 i.e., normal_phase
    scenes["s1/plots/p1/min_value"] = (config.initialCondition == 2 || config.initialCondition == 0)
                                       ? 0.0 :  matep.gap_A_td(config.p, config.T) * (1. + config.clamp_bias_gapMin);
    scenes["s1/plots/p1/max_value"] = matep.gap_B_td(config.p, config.T) * (1. + config.clamp_bias_gapMax);
    
    //scenes["s1/image_prefix"] = "ascent_output_render_gapA";
    scenes["s1/renders/r1/image_prefix"] = "gapA_t-%04d";
    scenes["s1/renders/r1/camera/azimuth"] = config.camera_azi/*35.0*/;
    scenes["s1/renders/r1/camera/elevation"] = config.camera_ele/*30.0*/;

    /* >>>>>>>>>>>>>>   2nd scene    <<<<<<<<<<<<< */
    
    // our 2nd scene (named 's5') will render the field 'feDensityOrdered'
    // to the file out_scene_ex1_render_var1.png
    scenes["s5/plots/p1/type"] = "pseudocolor";
    scenes["s5/plots/p1/field"] = "feDensityOrdered";

    // color map clamping. min_value will be set to 0.0 if initialCondtion is 2 i.e., normal_phase
    scenes["s5/plots/p1/min_value"] = 0.0 /*(config.initialCondition == 2 || config.initialCondition == 0)
					     ? 0.0 :  matep.gap_A_td(config.p, config.T) * (1. + config.clamp_bias_gapMin)*/;
    scenes["s5/plots/p1/max_value"] = 5.0;/*matep.gap_B_td(config.p, config.T) * (1. + config.clamp_bias_gapMax)*/;
    
    scenes["s5/renders/r1/image_prefix"] = "feDensity_t-%04d";
    scenes["s5/renders/r1/camera/azimuth"] = config.camera_azi/*35.0*/;
    scenes["s5/renders/r1/camera/elevation"] = config.camera_ele/*30.0*/;
    

    /* >>>>>>>>>>>>>> pipleline clip <<<<<<<<<<<<< */
    
    conduit::Node &add_act2 = actions.append();
    add_act2["action"] = "add_pipelines";
    conduit::Node &pipelines = add_act2["pipelines"];
    
    // pipeline 1
    pipelines["pl1/f1/type"] = "clip";
    conduit::Node &clip_params = pipelines["pl1/f1/params"];
    //clip_params["field"] = "gapAOrdered";
    clip_params["topology"] = "topo";
    clip_params["plane/point/x"] = 40.;
    clip_params["plane/point/y"] = 32.;
    clip_params["plane/point/z"] = 32.;
    clip_params["plane/normal/x"] = 0.;
    clip_params["plane/normal/y"] = 0.;
    clip_params["plane/normal/z"] = -1.;

    /*--------------------------------------------------*/
    /*conduit::Node &add_act3 = actions.append();
    add_act3["action"] = "add_extracts";

    conduit::Node &extracts = add_act3["extracts"];
    extracts["e1/type"] = "relay";
    extracts["e1/pipeline"]  = "pl1";
    extracts["e1/params/path"] = "gapA_clip";
    extracts["e1/params/protocol"] = "blueprint/mesh/hdf5";*/
    /*--------------------------------------------------*/

    scenes["s2/plots/p1/type"] = "pseudocolor";
    scenes["s2/plots/p1/pipeline"] = "pl1";
    scenes["s2/plots/p1/field"] = "gapAOrdered";

    // color map clamping. min_value will be set to 0.0 if initialCondtion is 2 i.e., normal_phase
    scenes["s2/plots/p1/min_value"] = (config.initialCondition == 2 || config.initialCondition == 0)
                                       ? 0.0 :  matep.gap_A_td(config.p, config.T) * (1. + config.clamp_bias_gapMin);
    scenes["s2/plots/p1/max_value"] = matep.gap_B_td(config.p, config.T) + config.clamp_bias_gapMax;
    
    scenes["s2/renders/r1/image_prefix"] = "gapA-clip_t-%04d";
    scenes["s2/renders/r1/camera/azimuth"] = config.camera_azi/*35.0*/;
    scenes["s2/renders/r1/camera/elevation"] = config.camera_ele/*30.0*/;

    /* >>>>>>>>>>> pipleline isosurfece <<<<<<<<<<<<< */

    /*Node &add_act3 = actions.append();
    add_act3["action"] = "add_pipelines";
    Node &pipelines2 = add_act3["pipelines"];
    // create a  pipeline (pl1) with a contour filter (f1)
    pipelines2["pl1/f1/type"] = "contour";*/
    pipelines["pl2/f1/type"] = "contour";

    // extract contours where gapAOrdered equals 2.0 to 3.0
    conduit::Node &contour_params = pipelines["pl2/f1/params"];
    contour_params["field"] = "gapAOrdered";
    //double iso_vals[21] = {3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4.0, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 4.5, 4.55, 4.6, 4.65, 4.7};
    double iso_vals[2] = {3.45, 3.55};
    //double iso_vals[11] = {2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0};
    //double iso_vals[12] = {2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9};

    contour_params["iso_values"].set(iso_vals,2);

    scenes["s3/plots/p1/type"] = "pseudocolor";
    scenes["s3/plots/p1/pipeline"] = "pl2";
    scenes["s3/plots/p1/field"] = "gapAOrdered";
    scenes["s3/renders/r1/image_prefix"] = "gapA-iso_t-%04d";

    double box_bounds[6] = {0.0, config.lx * config.dx, 0.0, config.ly * config.dx, 0.0, config.lz * config.dx};
    scenes["s3/renders/r1/dataset_bounds"].set(box_bounds,6);
    scenes["s3/renders/r1/render_bg"] = "true";
    
    scenes["s3/renders/r1/camera/azimuth"] = config.camera_azi/*35.0*/;
    scenes["s3/renders/r1/camera/elevation"] = config.camera_ele/*30.0*/;

    /* >>>>>>>>>>>>> extract hdf5 <<<<<<<<<<<<<< */
    
    conduit::Node &add_act3 = actions.append();
    add_act3["action"] = "add_extracts";

    conduit::Node &extracts = add_act3["extracts"];
    extracts["e1/type"] = "relay";
    extracts["e1/params/path"] = "sim-data";
    extracts["e1/params/protocol"] = "blueprint/mesh/hdf5";

    extracts["e1/params/fields"].append().set("gapAOrdered");
    //extracts["e1/params/fields"].append("gapAOrdered");
    extracts["e1/params/fields"].append().set("feDensityOrdered");
    
    extracts["e1/params/fields"].append().set("u11Ordered");
    extracts["e1/params/fields"].append().set("u12Ordered");
    extracts["e1/params/fields"].append().set("u13Ordered");
    extracts["e1/params/fields"].append().set("u21Ordered");
    extracts["e1/params/fields"].append().set("u22Ordered");
    extracts["e1/params/fields"].append().set("u23Ordered");
    extracts["e1/params/fields"].append().set("u31Ordered");
    extracts["e1/params/fields"].append().set("u32Ordered");
    extracts["e1/params/fields"].append().set("u33Ordered");

    extracts["e1/params/fields"].append().set("v11Ordered");
    extracts["e1/params/fields"].append().set("v12Ordered");
    extracts["e1/params/fields"].append().set("v13Ordered");
    extracts["e1/params/fields"].append().set("v21Ordered");
    extracts["e1/params/fields"].append().set("v22Ordered");
    extracts["e1/params/fields"].append().set("v23Ordered");
    extracts["e1/params/fields"].append().set("v31Ordered");
    extracts["e1/params/fields"].append().set("v32Ordered");
    extracts["e1/params/fields"].append().set("v33Ordered");
   
    /* >>>>>>>>>>>>>>   3slice  <<<<<<<<<<<<<<<< */
    pipelines["pl3/f1/type"] = "3slice";

    conduit::Node &slice3_params = pipelines["pl3/f1/params"];
    slice3_params["x_offset"] = -0.531f;
    slice3_params["y_offset"] = -0.6875f;
    slice3_params["z_offset"] = -0.21f;

    scenes["s4/plots/p1/type"] = "pseudocolor";
    scenes["s4/plots/p1/pipeline"] = "pl3";
    scenes["s4/plots/p1/field"] = "gapAOrdered";
    scenes["s4/renders/r1/image_prefix"] = "gapA-3slice_t-%04d";

    scenes["s4/plots/p1/min_value"] = 0.;
    scenes["s4/plots/p1/max_value"] = matep.gap_B_td(config.p, config.T) + config.clamp_bias_gapMax;
  
    scenes["s4/renders/r1/camera/azimuth"] = config.camera_azi/*35.0*/;
    scenes["s4/renders/r1/camera/elevation"] = config.camera_ele/*30.0*/;
  

    // print our full actions tree
    hila::out0 << actions.to_yaml() << '\n' << std::endl;
}

void he3sim::insitu_hdf5xdmf(){

  const std::string fname = "rank_xmls/" + config.xmf2_fname + "_" + std::to_string(hila::myrank()) + ".xml";
  //config.xdmf_out.open(config.xmf2_fname, std::ios::out);
  config.xdmf_out.open(fname, std::ios::out);
  const long dim_0 = lattice.mynode.size[0] + 2,
             dim_1 = lattice.mynode.size[1] + 2,
             dim_2 = lattice.mynode.size[2] + 2;

  unsigned int n;
  
    config.xdmf_out   << "<Grid Name=\"sim-data\" Type=\"Uniform\">\n"
                      << "  <Topology name=\"topo\" TopologyType=\"3DRectMesh\" Dimensions=\""
		      << dim_2 << " " << dim_1 << " " << dim_0 << "\"" << ">" << "\n"
                      << "  </Topology>\n"
                      << "  <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n"
                      << "    <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">"
		      << "\n"
                      << "      "
                      << ((lattice.mynode.min[0] - 1) * config.dx) << " "
                      << ((lattice.mynode.min[1] - 1) * config.dx) << " "
                      << ((lattice.mynode.min[2] - 1) * config.dx) << "\n"
                      << "    </DataItem>\n"
                      << "    <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n"
                      << "      "
                      << config.dx << " " << config.dx << " " << config.dx << "\n"
                      << "    </DataItem>\n"
                      << "  </Geometry>\n"
                      << "  <Attribute Name=\"gapA\" AttributeType=\"Scalar\" Center=\"Node\">\n"
                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	              << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank() << ".hdf5:/fields/gapAOrdered/values"
		      << "\n"
                      << "   </DataItem>\n"
                      << "  </Attribute>\n"
                      << "  <Attribute Name=\"feDensity\" AttributeType=\"Scalar\" Center=\"Node\">\n"
                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	              << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank() << ".hdf5:/fields/feDensityOrdered/values"
		      << "\n"
                      << "   </DataItem>\n"
                      << "  </Attribute>"
                      << "\n";

    for (n = 0; n<=8; ++n){
    config.xdmf_out   << "  <Attribute Name=\"u"
	              << std::to_string(n/3u + 1)
		      << std::to_string(n%3u + 1) << "\""
		      << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	              << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		      << ".hdf5:/fields/u"
		      << std::to_string(n/3u + 1)
		      << std::to_string(n%3u + 1) << "Ordered/values"
		      << "\n"
                      << "   </DataItem>\n"
                      << "  </Attribute>"
	              << "\n"
	              << "  <Attribute Name=\"v"
	              << std::to_string(n/3u + 1)
		      << std::to_string(n%3u + 1) << "\""
		      << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	              << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
		      << ".hdf5:/fields/v"
		      << std::to_string(n/3u + 1)
		      << std::to_string(n%3u + 1) << "Ordered/values"
		      << "\n"
                      << "   </DataItem>\n"
                      << "  </Attribute>"
	              << "\n";
    }

     
    config.xdmf_out   << "  <Attribute Name=\"vtkGhostType\" AttributeType=\"Scalar\" Center=\"Cell\">\n"
                      << "   <DataItem Format=\"HDF\" DataType=\"UChar\" Dimensions=\""
                      << dim_0 - 1 << " " << dim_1 - 1 << " " << dim_2 - 1 << "\"" << ">" << "\n"
 	              << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank() << ".hdf5:/fields/ascent_ghosts/values"
		      << "\n"   
                      << "   </DataItem>\n"
                      << "  </Attribute>\n"
                      << "</Grid>"
                      << "\n";

  config.xdmf_out.close();

}

void he3sim::insitu_execute() {

  //gapA[ALL] = /*phi[X].abs();*/ // A.degger() A. Trace() stuff......
    gapA[ALL] = real(sqrt((A[X]*A[X].dagger()).trace()));  

    /*------------------------------------------------------------*/
    /*--------------------     feDensity      --------------------*/
    real_t ebfe=fmin(matep.f_A_td(config.p, config.T), matep.f_B_td(config.p, config.T));
    feDensity[ALL] = 0;
    onsites(ALL) {
      Complex<double> a(0),b2(0),b3(0),b4(0),b5(0);
      //Complex<double> kin(0);
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
      //kin = (pi[X]*pi[X].dagger()).trace();
      
      foralldir(j) foralldir (k) foralldir(al){
	k1 += (A[X + k].column(j) - A[X - k].column(j)).e(al)
	      * (A[X + k].conj().column(j) - A[X - k].conj().column(j)).e(al)/(4.0*config.dx*config.dx);
	k2 += (A[X + j].column(j) - A[X - j].column(j)).e(al)
	      * (A[X + k].conj().column(k) - A[X - k].conj().column(k)).e(al)/(4.0*config.dx*config.dx);
	k3 += (A[X + k].column(j) - A[X - k].column(j)).e(al)
	      * (A[X + j].conj().column(k) - A[X - j].conj().column(k)).e(al)/(4.0*config.dx*config.dx);
      }
    // question here : how about imagnary part of k1 + k2 + k3 +bfe
      feDensity[X] = real(k1 + k2 + k3 + bfe);
    } //onsite(All) end here
    
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /*----------------     A matrix elements ---------------------*/
    u11[ALL] = A[X].e(0,0).re; v11[ALL] = A[X].e(0,0).im;
    u12[ALL] = A[X].e(0,1).re; v12[ALL] = A[X].e(0,1).im;
    u13[ALL] = A[X].e(0,2).re; v13[ALL] = A[X].e(0,2).im;
    u21[ALL] = A[X].e(1,0).re; v21[ALL] = A[X].e(1,0).im;
    u22[ALL] = A[X].e(1,1).re; v22[ALL] = A[X].e(1,1).im;
    u23[ALL] = A[X].e(1,2).re; v23[ALL] = A[X].e(1,2).im;
    u31[ALL] = A[X].e(2,0).re; v31[ALL] = A[X].e(2,0).im;
    u32[ALL] = A[X].e(2,1).re; v32[ALL] = A[X].e(2,1).im;
    u33[ALL] = A[X].e(2,2).re; v33[ALL] = A[X].e(2,2).im;       

    /*------------------------------------------------------------*/
        
  gapA.copy_local_data_with_halo(gapAOrdered);
  feDensity.copy_local_data_with_halo(feDensityOrdered);
  u11.copy_local_data_with_halo(u11Ordered); v11.copy_local_data_with_halo(v11Ordered);
  u12.copy_local_data_with_halo(u12Ordered); v12.copy_local_data_with_halo(v12Ordered);
  u13.copy_local_data_with_halo(u13Ordered); v13.copy_local_data_with_halo(v13Ordered);
  u21.copy_local_data_with_halo(u21Ordered); v21.copy_local_data_with_halo(v21Ordered);
  u22.copy_local_data_with_halo(u22Ordered); v22.copy_local_data_with_halo(v22Ordered);
  u23.copy_local_data_with_halo(u23Ordered); v23.copy_local_data_with_halo(v23Ordered);
  u31.copy_local_data_with_halo(u31Ordered); v31.copy_local_data_with_halo(v31Ordered);
  u32.copy_local_data_with_halo(u32Ordered); v32.copy_local_data_with_halo(v32Ordered);
  u33.copy_local_data_with_halo(u33Ordered); v33.copy_local_data_with_halo(v33Ordered);       

  /* ToDo list :
   * > orbital vectors;
   * > spin vectors;
   * > mass current;
   * > spin currents;
   * x free energy density (done)
   * > ...
   */
    
  insitu.execute(actions);
}

void he3sim::insitu_initialize() {

    latticeVolumeWithGhost =
        (lattice.mynode.size[0] + 2) * (lattice.mynode.size[1] + 2) * (lattice.mynode.size[2] + 2);
    
    latticeVolume =
      (lattice.mynode.size[0]) * (lattice.mynode.size[1]) * (lattice.mynode.size[2]);

    gapAOrdered.reserve(latticeVolumeWithGhost);
    feDensityOrdered.reserve(latticeVolumeWithGhost);

    u11Ordered.reserve(latticeVolumeWithGhost); v11Ordered.reserve(latticeVolumeWithGhost);
    u12Ordered.reserve(latticeVolumeWithGhost); v12Ordered.reserve(latticeVolumeWithGhost);
    u13Ordered.reserve(latticeVolumeWithGhost); v13Ordered.reserve(latticeVolumeWithGhost);
    u21Ordered.reserve(latticeVolumeWithGhost); v21Ordered.reserve(latticeVolumeWithGhost);
    u22Ordered.reserve(latticeVolumeWithGhost); v22Ordered.reserve(latticeVolumeWithGhost);
    u23Ordered.reserve(latticeVolumeWithGhost); v23Ordered.reserve(latticeVolumeWithGhost);
    u31Ordered.reserve(latticeVolumeWithGhost); v31Ordered.reserve(latticeVolumeWithGhost);
    u32Ordered.reserve(latticeVolumeWithGhost); v32Ordered.reserve(latticeVolumeWithGhost);
    u33Ordered.reserve(latticeVolumeWithGhost); v33Ordered.reserve(latticeVolumeWithGhost);    
    
    // One more point in each direction, but cell data (Npts - 1 cells)
    auto ghostNX = lattice.mynode.size[0] + 2 - 1;
    auto ghostNY = lattice.mynode.size[1] + 2 - 1;
    auto ghostNZ = lattice.mynode.size[2] + 2 - 1;

    ghostVolume = ghostNX * ghostNY * ghostNZ;
    ghostCellsMask = (unsigned char *)memalloc(ghostVolume * sizeof(unsigned char));

    long long counter = 0;
    unsigned char Mask = 0;
    for (auto k = 0; k < ghostNZ; k++) {
        for (auto j = 0; j < ghostNY; j++) {
            for (auto i = 0; i < ghostNX; i++) {
                bool kGhostFlag = (k == 0);
                bool jGhostFlag = (j == 0);
                bool iGhostFlag = (i == 0);
                Mask = (iGhostFlag || jGhostFlag || kGhostFlag);
                ghostCellsMask[counter] = Mask;
                counter++;
            }
        }
    }

    insitu_createMesh();

    ascent_options["mpi_comm"] = MPI_Comm_c2f(lattice.mpi_comm_lat);
    ascent_options["runtime/type"] = "ascent";
#if defined CUDA
    ascent_options["runtime/vtkm/backend"] = "cuda";
    ascent_options["cuda/init"] = "false";
#endif
    ascent_options["timings"] = "false";

    // *************************
    insitu.open(ascent_options);
    insitu.publish(mesh);
    insitu_defineActions();
}

void he3sim::insitu_close() {
    insitu.close();
}
#endif

/*--------------------------------------------------------------------------------*/
/* In situ functions end at here */
/*--------------------------------------------------------------------------------*/


void he3sim::next() {

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
    // deltaPi[X] += (1.0/(4.0*config.dx*config.dx)) * (A[X + e_x + e_x] + A[X - e_x - e_x]
    // 						     + A[X + e_y + e_y] + A[X - e_y - e_y]
    // 						     + A[X + e_z + e_z] + A[X - e_z - e_z]
    // 						     - 6.0*A[X]);

    /*
     * debug gradient updates: 
     */
    deltaPi[X] +=  (1.0 / (config.dx * config.dx)) * (A[X + e_x] + A[X - e_x]
						      + A[X + e_y] + A[X - e_y]
						      + A[X + e_z] + A[X - e_z]
						      - 6.0 * A[X]);
  }

    //onsites (ALL) {deltaPi[X] *= config.dt;} // I think that this is the problem, multiplication with respect to dt

  if (t < config.tdif)
    {
      pi[ALL] = deltaPi[X]/(config.difFac);
      t += config.dt/config.difFac;
    }
  else if (t < config.tdis && config.gamma > 0 )
    {
      pi[ALL] = pi[X] + (deltaPi[X] - 2.0 * config.gamma * pi[X])*config.dt;
      t += config.dt;
    }
  else
    {
      pi[ALL] = pi[X] + deltaPi[X]*config.dt;
      t += config.dt;
    }
  

  next_timer.stop();

}

void he3sim::next_bath() {

  Matep MP;
  real_t Tp[2];
  update_Tp(t, Tp);

  static hila::timer next_timer("timestep");
  Field<phi_t> deltaPi;
  Field<Vector<3,Complex<real_t>>> djAaj;

  real_t gapa = MP.gap_A_td(Tp[1], Tp[0]);
  real_t gapb = MP.gap_B_td(Tp[1], Tp[0]);

  real_t tb = Tp[0]/tc;
  real_t sig = sqrt(2.0*tb*config.gamma);
  real_t ep2 = 1.0-exp(-2.0*config.gamma*config.dt) ;

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
  else if (t < config.tdis && config.gamma > 0 )
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

}

int main(int argc, char **argv) {
    he3sim sim;
    const std::string output_fname = sim.allocate("sim_params.txt", argc, argv);
    sim.update_params();
    sim.initialize();
    int stepspos;
    
    int steps =
        (sim.config.tEnd - sim.config.tStats) /
        (sim.config.dt * sim.config.nOutputs); // number of steps between printing stats
    if (steps == 0)
        steps = 1;

    /*if (sim.config.positions == 1)
      {
	stepspos =  (sim.config.tEnd - sim.config.tStats) /
	  (sim.config.dt * sim.config.npositionout);
	if (stepspos == 0)
        stepspos = 1;
	}*/
    /*if(sim.config.rank_stream_switch == 1) 
      {
       const std::string rankOut_fname = "rank-Output/AGap_rank_" + std::to_string(hila::myrank()) + ".dat";
       sim.config.rank_Output.open(rankOut_fname, std::ios::out);
       }*/
	
    int stat_counter = 0;

    if (hila::myrank() == 0) {
        sim.config.stream.open(output_fname, std::ios::out);
    }
    
    sim.insitu_hdf5xdmf();
    if (hila::myrank() == 0) {sim.write_xdmf();}
    
#if defined USE_ASCENT
    hila::out0 << "using ascent" << "\n\n\n";
    sim.insitu_initialize();
#endif    
      
    // on gpu the simulation timer is fake, because there's no sync here.  
    // BUt we want to avoid unnecessary sync anyway.
    static hila::timer run_timer("Simulation time"), meas_timer("Measurements");
    run_timer.start();

    // lattice coordinates output
    sim.latticeCoordinate_output();
    sim.insitu_hdf5xdmf();
    // output steps just before while loop
    // std::cout << "\n" << "steps is " << steps << " with nOutput " << sim.config.nOutputs << "\n";
    
    //auto tildephi = sim.phi;
    while (sim.t < sim.config.tEnd) {
        if (sim.t >= sim.config.tStats) {
	  
            if (stat_counter % steps == 0) {

              if(sim.config.rank_stream_switch == 1) 
               {
		 const std::string rankOut_fname = "rank-Output/AGap_rank_" + std::to_string(hila::myrank()) + "_t" + std::to_string(sim.t) + ".dat";
                 sim.config.rank_Output.open(rankOut_fname, std::ios::out);
               }            	      
	      meas_timer.start();
	      hila::out0 << "Writing output at time " << sim.t 
			 << ", stat_counter is " << stat_counter
                         << ", t point is " << sim.t/sim.config.dt
			 << ", steps is " << steps
			 << "\n";
	      sim.write_moduli();
	      //sim.write_energies();
	      //sim.write_phases();
	      //sim.write_gapA_rankOut();

#if defined USE_ASCENT
              sim.insitu_execute();
#endif
	      
	      meas_timer.stop();

	      /* close rank-output :*/
              sim.config.rank_Output.close();
            }
	    if (sim.config.positions == 1)
	      {
		if (stat_counter % stepspos == 0)
		  {
		  //sim.write_A_matrix_positions();
                   sim.write_positions();
		  }
	      }
            stat_counter++;
        }
	sim.update_params();
	if (sim.config.useTbath == 1)
	  {
	    sim.next_bath();}
	else
	  {
	    sim.next();}
    }
    run_timer.stop();

#if defined USE_ASCENT
    hila::out0 << "ascent closed" << "\n\n\n";
    sim.insitu_close();
#endif    

    if (hila::myrank() == 0) {
        sim.config.stream.close();
    }

    //    sim.config.rank_Output.close();


    hila::finishrun();
    return 0;
}

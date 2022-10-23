#define _USE_MATH_DEFINES
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
//#include <math.h>
#include <assert.h>


#include "plumbing/hila.h"
#include "plumbing/fft.h"

#include "matep.hpp"

// Definition of the fiedl that we will use
using real_t = float;   // or double
using phi_t = Matrix<3,3,Complex<real_t>>;     // saves the trouble of writing this every time



/// Container for simulation parameters and methods
class scaling_sim{

public:
  scaling_sim() = default;
  const std::string allocate(const std::string &fname, int argc, char **argv);
  void initialize();
  void update_params();
  void update_Tp (real_t t, real_t Tp[2]);
  void write_moduli();
  void write_energies();
  void write_positions();
  void next();
  int  boundary_point(CoordinateVector e);
  
  Field<phi_t> A;
  Field<phi_t> pi;

  real_t t;

  real_t tc = 0;

  std::vector<real_t> t_v;
  std::vector<real_t> T_v;
  std::vector<real_t> p_v;
  
    struct config {
      int l;
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
      int seed;
      real_t IniMod;
      real_t Inilc;
      
      int item;
      real_t T;
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
      int positions;
      std::fstream stream_positions;
      int boundary_conditions;
    } config;
};

const std::string scaling_sim::allocate(const std::string &fname, int argc,
                                        char **argv) {

   hila::initialize(argc, argv);

    hila::input parameters(fname);
    config.l = parameters.get("N");
    config.dx = parameters.get("dx");
    config.dtdxRatio = parameters.get("dtdxRatio");
    config.tStart = parameters.get("tStart");
    config.tEnd = parameters.get("tEnd");
    config.tdif = parameters.get("tdif");
    config.difFac = parameters.get("difFac");
    config.tdis = parameters.get("tdis");
    config.gamma = parameters.get("gamma");
    config.initialCondition = parameters.get("initialCondition");
    config.seed = parameters.get("seed");
    config.IniMod = parameters.get("IniMod");
    config.Inilc = parameters.get("Inilc");
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
      config.T = parameters.get("T");
      config.p = parameters.get("p");
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
    config.tStats = parameters.get("tStats");
    config.nOutputs = parameters.get("nOutputs");

    const std::string output_file = parameters.get("output_file");
    
    config.dt = config.dx * config.dtdxRatio;
    t = config.tStart;

    CoordinateVector box_dimensions = {config.l, config.l, config.l};
    lattice.setup(box_dimensions);
    hila::seed_random(config.seed);

    config.positions = 0;
    config.boundary_conditions=0;

    return output_file;
}


void scaling_sim::initialize() {

  Matep MP;
  real_t Tp[2];
  update_Tp(t, Tp);
  
  int N = config.l;
  real_t dx = config.dx;

  hila::out0<<"R: "<<hila::random()<<"\n";
  
  switch (config.initialCondition) {
    
  case 1: {
    pi = 0;
    real_t gap = MP.gap_B_td(Tp[1], Tp[0]);
    onsites(ALL) {
        A[X].gaussian_random();
      A[X] = gap * A[X]/A[X].norm();
    }

    hila::out0 << "Components randomly created \n";

    break;
    }
  case 2: {
    auto kA = A;
    real_t gap = config.IniMod; //MP.gap_B_td(Tp[1], Tp[0]);
    real_t lc = config.Inilc;//1.0/sqrt(abs(config.alpha));

    hila::out0 << "Correlation length in ICs: "<< lc <<"\n";
    
    onsites (ALL) {
            real_t constant = pow(gap, 2.0) * pow(2.0 * M_PI, 1.5) * pow(lc, 3.0)/(9.0 * N * N * N * dx * dx * dx);
            real_t kSqu;
            real_t std;
            kSqu = 0.0;
            auto k = X.coordinates();

            foralldir (d) { kSqu += pow(sin(M_PI * k.e(d) / N), 2.0); }
            kSqu *= pow(2.0 / dx, 2.0);

            if (kSqu > 0.0) {
                std = sqrt(0.5 * constant *
                           exp(-0.5 * kSqu * lc * lc));
		        kA[X].gaussian_random(std);

            } else {
	      kA[X]=0;
	      /*foralldir(d1) foralldir(d2) {
                  kA[X].e(d1,d2).re=0.0;
                  kA[X].e(d1,d2).im=0.0;
		  }*/
            }
    }

        FFT_field(kA, A, fft_direction::back);

	real_t rre=hila::gaussrand();
	real_t rim=hila::gaussrand();
	
	onsites (ALL){
	  foralldir(d1) foralldir(d2) {
	    A[X].e(d1,d2).re=1.0;
	    A[X].e(d1,d2).im=1.0;
	    }
	  A[X]=gap*A[X]/A[X].norm();
	}
	
        pi[ALL] = 0;

        hila::out0 << "k space generation \n";

        break;
  }
  case 3: {
    pi = 0;
    real_t gap = MP.gap_B_td(Tp[1], Tp[0]);
    hila::out0<<"Gap B: "<<gap<<"\n";
    onsites(ALL) {
      A[X] = gap/sqrt(3.0);   // there was 1/A[X].norm(), but this would have been 1/sqrt(3)
    }

    hila::out0 << "Pure B phase \n";

    break;
    }
  case 4: {
    pi = 0;
    real_t gap = MP.gap_A_td(Tp[1], Tp[0]);
    hila::out0<<"Gap A: "<<gap<<"\n";
    onsites(ALL) {
        A[X] = 0;
        A[X].e(0,2) = 1;
        A[X].e(1,2) = I;
      
      A[X] = gap * A[X]/sqrt(2.0);
    }

    hila::out0 << "Pure A phase \n";

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

void scaling_sim::update_Tp (real_t t, real_t Tp[2]) {

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
void scaling_sim::update_params() {

  Matep MP;
  real_t Tp[2];

  
  if (config.item ==1 && t == config.tStart){
    tc = MP.Tcp_mK(config.p);
    config.alpha = MP.alpha_td(config.p, config.T);
    config.beta1 = MP.beta1_td(config.p, config.T);
    config.beta2 = MP.beta2_td(config.p, config.T);
    config.beta3 = MP.beta3_td(config.p, config.T);
    config.beta4 = MP.beta4_td(config.p, config.T);
    config.beta5 = MP.beta5_td(config.p, config.T);
    hila::out0 << config.alpha << " " << config.beta1 << " " << config.beta2 << " " << config.beta3 << " " << config.beta4 << " " << config.beta5;
  }
  else if (config.item ==2){
    update_Tp(t, Tp);
    tc = MP.Tcp_mK(Tp[1]);
    config.alpha = MP.alpha_td(Tp[1], Tp[0]);
    config.beta1 = MP.beta1_td(Tp[1], Tp[0]);
    config.beta2 = MP.beta2_td(Tp[1], Tp[0]);
    config.beta3 = MP.beta3_td(Tp[1], Tp[0]);
    config.beta4 = MP.beta4_td(Tp[1], Tp[0]);
    config.beta5 = MP.beta5_td(Tp[1], Tp[0]);
  }
}

void scaling_sim::write_moduli() {

   // real_t a = scaleFactor(t);

    double Amod = 0.0;
    double pimod = 0.0;

    real_t Tp[2];

    hila::set_allreduce(false);
    onsites (ALL) {
        Amod += A[X].norm();
        pimod += pi[X].norm();
    }

        
    update_Tp(t, Tp);
    
    if (hila::myrank() == 0) {
        config.stream << t << " "
	              << tc << " "
	              << Tp[0] << " " << Tp[1] << " "
		      << config.alpha << " " << config.beta1 << " " <<	config.beta2 << " " <<	config.beta3 << " " <<	config.beta4 << " " <<	config.beta5 << " "
                      << Amod / lattice.volume() << " " << pimod / lattice.volume()
                      << " ";
    }
}

void scaling_sim::write_energies() {

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
	k1 += (A[X + k].column(j) - A[X - k].column(j)).e(al) * (A[X + k].conj().column(j) - A[X - k].conj().column(j)).e(al)/(4.0*config.dx*config.dx);
	k2 += (A[X + j].column(j) - A[X - j].column(j)).e(al) * (A[X + k].conj().column(k) - A[X - k].conj().column(k)).e(al)/(4.0*config.dx*config.dx);
	k3 += (A[X + k].column(j) - A[X - k].column(j)).e(al) * (A[X + j].conj().column(k) - A[X - j].conj().column(k)).e(al)/(4.0*config.dx*config.dx);
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
	              << sumb5_we.re / vol << " " << sumb5_we.im / vol << "\n";
    }

       hila::out0<< "Energy done \n";

    //   double sumar = 0.0;
    //   double sumai = 0.0;
    //   double sumb1 = 0.0;
    //   double sumb2r = 0.0;
    //   double sumb2i = 0.0;
    //   double sumb3r = 0.0;
    //   double sumb3i = 0.0;
    //   double sumb4r = 0.0;
    //   double sumb4i = 0.0;
    //   double sumb5r = 0.0;
    //   double sumb5i = 0.0;

    //   hila::set_allreduce(false);
    //   onsites (ALL) {

    //     sumar += config.alpha * ((A[X]*A[X].dagger()).trace()).re;
    //     sumai += config.alpha * ((A[X]*A[X].dagger()).trace()).im;
        
    //     sumb1 += config.beta1 * ((A[X]*A[X].transpose()).trace()).squarenorm();

    //     sumb2r += config.beta2 * ((A[X]*A[X].dagger()).trace()*(A[X]*A[X].dagger()).trace()).re;
    //     sumb2i += config.beta2 * ((A[X]*A[X].dagger()).trace()*(A[X]*A[X].dagger()).trace()).im;

    //     sumb3r += config.beta3 * ((A[X]*A[X].transpose()*A[X].conj()*A[X].dagger()).trace()).re;
    //     sumb3i += config.beta3 * ((A[X]*A[X].transpose()*A[X].conj()*A[X].dagger()).trace()).im;

    //     sumb4r += config.beta4 * ((A[X]*A[X].dagger()*A[X]*A[X].dagger()).trace()).re;
    //     sumb4i += config.beta4 * ((A[X]*A[X].dagger()*A[X]*A[X].dagger()).trace()).im;

    //     sumb5r += config.beta5 * ((A[X]*A[X].dagger()*A[X].conj()*A[X].transpose()).trace()).re;
    //     sumb5i += config.beta5 * ((A[X]*A[X].dagger()*A[X].conj()*A[X].transpose()).trace()).im;
    //   }

    //   if (hila::myrank() == 0) {
    //     double vol = (double)config.l * config.l * config.l;
    //     config.stream << sumar / vol << " "<< sumai / vol << " "
    // 		  << sumb1 / vol << " "
    // 		  << sumb2r / vol << " " << sumb2i / vol << " "
    // 		  << sumb3r / vol << " " << sumb3i / vol << " "
    //  		  << sumb4r / vol << " " << sumb4i / vol << " "
    // 		  << sumb5r / vol << " " << sumb5i / vol << "\n";
    //   }
}

void scaling_sim::write_positions() {

  Matep MP;

      real_t Tp[2];

      update_Tp(t, Tp);

      real_t ebfe=fmin(MP.f_A_td(Tp[1], Tp[0]),MP.f_B_td(Tp[1], Tp[0]));

                 
      if (hila::myrank() == 0)
        {
          config.stream_positions<< t << " " << MP.f_A_td(Tp[1], Tp[0])<<" "<<MP.f_B_td(Tp[1], Tp[0]) <<"\n";
	}
      
      hila::set_allreduce(false);
      onsites(ALL) {

	Complex<double> a(0),b2(0),b3(0),b4(0),b5(0);
	Complex<double> bfe(0);
	double b1 = 0;

	a = config.alpha * (A[X]*A[X].dagger()).trace();

	b1 = config.beta1 * ((A[X]*A[X].transpose()).trace()).squarenorm();

	b2 = config.beta2 * ((A[X]*A[X].dagger()).trace()*(A[X]*A[X].dagger()).trace());

	b3 = config.beta3 * ((A[X]*A[X].transpose()*A[X].conj()*A[X].dagger()).trace());

	b4 = config.beta4 * ((A[X]*A[X].dagger()*A[X]*A[X].dagger()).trace());

	b5 = config.beta5 * ((A[X]*A[X].dagger()*A[X].conj()*A[X].transpose()).trace());

	bfe = a + b1 + b2 + b3 + b4 + b5;

	if (bfe.re<-0.07)
	  {
	    config.stream_positions<< bfe <<"\n";
	    
	    }

	
	/*if ( hila::myrank() == 0)
	  {config.stream_positions<< bfe.re << " " << X.coordinate(e_x) << " " << X.coordinate(e_y) << " "<< X.coordinate(e_z) << " "
	       << A[X].e(0,0).re << " " << A[X].e(0,0).re
				  << A[X].e(0,1).re << " " << A[X].e(0,1).re
				  << A[X].e(0,2).re << " " << A[X].e(0,2).re
				  << A[X].e(1,0).re << " " << A[X].e(1,0).re
				  << A[X].e(1,1).re << " " << A[X].e(1,1).re
				  << A[X].e(1,2).re << " " << A[X].e(1,2).re
				  << A[X].e(2,0).re << " " << A[X].e(2,0).re
				  << A[X].e(2,1).re << " " << A[X].e(2,1).re
				  << A[X].e(2,2).re << " " << A[X].e(2,2).re
				  <<  "\n";
	  }*/
      
      }

}

int  scaling_sim::boundary_point(CoordinateVector e) {

  if (e[e_x]==0 or e[e_x]==(config.l-1) or e[e_y]==0 or e[e_y]==(config.l-1) or e[e_z]==0 or e[e_z]==(config.l-1))
    {
      return 1;
    }
  else
    {
      return 0;
    }
}

void scaling_sim::next() {


    static hila::timer next_timer("timestep");
    Field<phi_t> deltaPi;
    Field<Vector<3,Complex<real_t>>> djAaj;
    
    next_timer.start();

    onsites (ALL) {

      if (config.boundary_conditions == 1)
	{
	  if (boundary_point(X.coordinates())==1)
	    {A[X] = 0.0;}
	  else {A[X] += config.dt * pi[X];}
	}
      else {A[X] += config.dt * pi[X];}


      auto AxAt = A[X]*A[X].transpose();
      auto AxAd = A[X]*(A[X].dagger());
 
	
      deltaPi[X] = - config.alpha*A[X] - 2.0*config.beta1*AxAt.trace()*A[X].conj() 
	- 2.0*config.beta2*AxAd.trace()*A[X]
	- 2.0*config.beta3*(AxAt*A[X].conj()) - 2.0*config.beta4*(AxAd*A[X]) 
	  //- 2.0*config.beta5*(AxAd.conj()*A[X]);
	- 2.0*config.beta5*(A[X].conj()*A[X].transpose()*A[X]);

    //     deltaPi[X] = - config.alpha*A[X] - 2.0*config.beta1*A[X]*(A[X]*A[X].transpose()).trace() - 2.0*config.beta2*A[X]*(A[X]*A[X].dagger()).trace()
	//   - 2.0*config.beta3*(A[X]*A[X].transpose()*A[X]) - 2.0*config.beta4*(A[X]*A[X].dagger()*A[X]) - 2.0*config.beta5*(A[X].conj()*A[X].transpose()*A[X]);
    }

    onsites(ALL) {
        djAaj[X] = 0;
        foralldir(j) {
	  auto Xp = X.coordinates()+Direction(j);
	  auto Xm = X.coordinates()-Direction(j);
	  if (config.boundary_conditions == 1)
	    {
	      if (boundary_point(Xp)==1 and boundary_point(Xm)==1)
		{
		  djAaj[X] = djAaj[X];
		}
	      else if (boundary_point(Xp)==1)
		{
		  djAaj[X] +=  - A[X - j].column(j);
		}
	      else if (boundary_point(Xm)==1)
		{
		  djAaj[X] += A[X + j].column(j);
		}
	      else {djAaj[X] += A[X + j].column(j) - A[X - j].column(j);}
	    }
	  else {djAaj[X] += A[X + j].column(j) - A[X - j].column(j);}
	    
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
      if (config.boundary_conditions == 1)
	{
	  foralldir(j) {
	    auto Xp = X.coordinates()+Direction(j);
	    auto Xm = X.coordinates()-Direction(j);
	    if (boundary_point(Xp)==1 and boundary_point(Xm)==1)
	      {
		deltaPi[X] = deltaPi[X]; 
	      }
	    else if (boundary_point(Xp)==1)
	      {
		deltaPi[X] += (1.0/(config.dx*config.dx)) * (A[X -e_x]);
	      }
	    else if (boundary_point(Xm)==1)
	      {
		deltaPi[X] += (1.0/(config.dx*config.dx)) * (A[X -e_x]);
	      }
	    else
	      {
		deltaPi[X] += (1.0/(config.dx*config.dx)) * (A[X + e_x] + A[X -e_x]);
	      }
	    if (boundary_point(X.coordinates())==1)
	      {
		deltaPi[X] = deltaPi[X];
	      }
	    else
	      {
		deltaPi[X] += (1.0/(config.dx*config.dx)) * ( - 6.0*A[X]);
	      }
	  }
	}
      else
	{
	  deltaPi[X] += (1.0/(4.0*config.dx*config.dx)) * (A[X + e_x + e_x] + A[X - e_x - e_x]
						     + A[X + e_y + e_y] + A[X - e_y - e_y]
						     + A[X + e_z + e_z] + A[X - e_z - e_z]
						   - 6.0*A[X]);
	}
    }

    //onsites (ALL) {deltaPi[X] *= config.dt;} // I think that this is the problem, multiplication with respect to dt

      if (t < config.tdif)
	{
	  pi[ALL] = deltaPi[X]/(config.difFac);
	  t += config.dt/config.difFac;
	}
      else if (t < config.tdis && config.gamma > 0 )
	{
	  pi[ALL] = pi[X] + (deltaPi[X]*config.dt - 2.0 * config.gamma * pi[X]);
	  t += config.dt;
	}
      else
	{
	  onsites (ALL){pi[X] = pi[X] + deltaPi[X]*config.dt;}
	  t += config.dt;
	}


    next_timer.stop();

}

int main(int argc, char **argv) {
    scaling_sim sim;
    const std::string output_fname = sim.allocate("sim_params.txt", argc, argv);
    sim.update_params();
    sim.initialize();
    
    int steps =
        (sim.config.tEnd - sim.config.tStats) /
        (sim.config.dt * sim.config.nOutputs); // number of steps between printing stats
    if (steps == 0)
        steps = 1;
        
    int stat_counter = 0;

    if (hila::myrank() == 0) {
        sim.config.stream.open(output_fname, std::ios::out);
    }

    if (sim.config.positions == 1 and hila::myrank() == 0)
      {
	const std::string positions_fname = "positions_"+std::to_string(hila::myrank())+".dat";
	sim.config.stream_positions.open(positions_fname, std::ios::out);
      }
    
    // on gpu the simulation timer is fake, because there's no sync here.  
    // BUt we want to avoid unnecessary sync anyway.
    static hila::timer run_timer("Simulation time"), meas_timer("Measurements");
    run_timer.start();
    
    //auto tildephi = sim.phi;
    while (sim.t < sim.config.tEnd) {
        if (sim.t >= sim.config.tStats) {
            if (stat_counter % steps == 0) {
	      meas_timer.start();
	      hila::out0 << "Writing output at time " << sim.t << "\n";
	      if (sim.config.positions == 1){sim.write_positions();}
	      sim.write_moduli();
	      sim.write_energies();
	      meas_timer.stop();
            }
            stat_counter++;
        }
	sim.update_params();
        sim.next();
    }
    run_timer.stop();

    if (hila::myrank() == 0) {
        sim.config.stream.close();
	sim.config.stream_positions.close();
    }

    hila::finishrun();
    return 0;
}

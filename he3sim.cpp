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

// Definition of the fieldthat we will use
using real_t = float;                          // or double ?
using phi_t = Matrix<3,3,Complex<real_t>>;     // saves the trouble of writing this every time



// Container for simulation parameters and methods
class scaling_sim{

public:
  scaling_sim() = default;
  const std::string allocate(const std::string &fname, int argc, char **argv);
  void initialize();
  void initializeT();
  void initializep();
  void point_params(real_t T, real_t p, real_t beta[6]);
  void write_moduli();
  void write_energies();
  void write_positions();
  void write_phases();
  void next();
  void next_bath();
  void nextT();


  void write_A_matrix_positions();             // output A-matrix after certain time interval
  void latticeCoordinate_output();             // lattice coordinates output with same sequence of A.write() 
  
  Field<phi_t> A;
  Field<phi_t> pi;

  Field<real_t> T;
  Field<real_t> p;
  
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
      int seed;
      real_t IniMod;
      real_t Inilc;
      
      int item;
      int initialConditionT;
      real_t IniT;
      real_t ampT;
      real_t sigT;
      int initialConditionp;
      real_t Inip;
      //real_t T;
      //real_t p;
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
      int npositionout;
      int boundary_conditions;
      int useTbath;
      int write_phases;
      int write_eigen;
      int evolveT;
      real_t diffT;
    } config;
};

const std::string scaling_sim::allocate(const std::string &fname, int argc, char **argv)
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
    config.initialCondition = parameters.get_item("initialCondition",{"gaussrand", "kgaussrand","Bphase","Aphase"});
    config.seed = parameters.get("seed");
    config.IniMod = parameters.get("IniMod");
    config.Inilc = parameters.get("Inilc");
    config.initialConditionT = parameters.get_item("initialConditionT",{"constant","sine","hotspot"});
    if(config.initialConditionT == 0)
      {
	config.IniT = parameters.get("IniT");
      }
    else if (config.initialConditionT == 1)
      {
	config.IniT = parameters.get("IniT");
	config.ampT = parameters.get("ampT");
      }
    else if (config.initialConditionT == 2)
      {
	config.IniT = parameters.get("IniT");
        config.ampT = parameters.get("ampT");
	config.sigT = parameters.get("sigT");
      }
    config.initialConditionp = parameters.get_item("initialConditionp",{"constant"});
    config.Inip	= parameters.get("Inip");
    config.tStats = parameters.get("tStats");
    config.nOutputs = parameters.get("nOutputs");

    // output_file is the saving path of output file, which offered in congigration file
    const std::string output_file = parameters.get("output_file");

    config.positions = parameters.get_item("out_points",{"no", "yes"});
    if(config.positions==1)
      {
	config.npositionout = parameters.get("npositionout");
	config.write_phases = parameters.get_item("write_phases",{"no","yes"});
	config.write_eigen = parameters.get_item("write_eigen",{"no","yes"});
      }
    config.boundary_conditions = parameters.get_item("boundary_conditions",{"periodic", "AB", "PairBreaking"});
    config.useTbath = parameters.get_item("useTbath",{"no","yes"});
    config.evolveT = parameters.get_item("evolveT",{"no","yes"});
    if(config.evolveT ==1)
      {
	config.diffT = parameters.get("diffT");
      }
    
    config.dt = config.dx * config.dtdxRatio;
    t = config.tStart;

    // setup the hila lattice geometry 
    CoordinateVector box_dimensions = {config.lx, config.ly, config.lz};
    lattice.setup(box_dimensions);
    hila::seed_random(config.seed);


    return output_file;
}


void scaling_sim::initialize() {

  Matep MP;
   
  int Nx = config.lx;
  int Ny = config.ly;
  int Nz = config.lz;
  
  real_t dx = config.dx;

  real_t ttc = MP.Tcp_mK(config.Inip);
  hila::out0 <<"T_AB: "<<MP.tAB_RWS(config.Inip)*ttc<<"\n";
  
  switch (config.initialCondition) {
    
  case 0: {
    pi = 0;                            
    real_t gap = MP.gap_B_td(config.Inip, config.IniT);
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
    real_t gap = MP.gap_B_td(config.Inip, config.IniT);
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
  case 3: {
    pi = 0;
    real_t gap = MP.gap_A_td(config.Inip, config.IniT);
    hila::out0<<"Gap A: "<<gap<<"\n";
    onsites(ALL) {

        A[X] = 0;
        A[X].e(0,2).re = 1.0;         // what indice of e means? Is this A-phase?
        A[X].e(1,2).im = 1.0;
      
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

void scaling_sim::initializeT() {

  switch (config.initialCondition) {

  case 0: {
    
    onsites(ALL) {
      T[X] = config.IniT;
    }
    break;
  }
   
  case 1: {

    onsites(ALL){
      
      auto xcoord = X.coordinate(e_x);

      T[X] = config.IniT + config.ampT*sin(2.0*M_PI*xcoord/config.lx);
    
    }
    break;
  }

  case 2: {

    onsites(ALL){
      auto xcoord = X.coordinate(e_x);
      auto ycoord = X.coordinate(e_y);
      auto zcoord = X.coordinate(e_z);

      real_t r = sqrt((xcoord-config.lx/2.0)*(xcoord-config.lx/2.0)+(ycoord-config.ly/2.0)*(ycoord-config.ly/2.0)+(zcoord-config.lz/2.0)*(zcoord-config.lz/2.0));

      real_t expr=exp(-0.5*(r/config.sigT)*(r/config.sigT));

      T[X] = config.IniT + config.ampT*expr;
    }
    break;
  } 
  }
}

void scaling_sim::initializep() {

  switch (config.initialCondition) {

  case 0: {
    onsites(ALL) {

      p[X] = config.Inip;

    }
  }
  } 
}

void scaling_sim::point_params(real_t T, real_t p, real_t beta[6]) {

  Matep MP;
  
  beta[0] = MP.alpha_td(p, T);
  beta[1] = MP.beta1_td(p, T);
  beta[2] = MP.beta2_td(p, T);
  beta[3] = MP.beta3_td(p, T);
  beta[4] = MP.beta4_td(p, T);
  beta[5] = MP.beta5_td(p, T);
}
  
void scaling_sim::write_moduli() {

   // real_t a = scaleFactor(t);

    double Amod = 0.0;
    double pimod = 0.0;

    //real_t Tp[2];

    hila::set_allreduce(false);
    onsites (ALL) {
        Amod += A[X].norm();
        pimod += pi[X].norm();
    }

        
    //update_Tp(t, Tp);
    
    if (hila::myrank() == 0) {
        config.stream << t << " "
                      << Amod / lattice.volume() << " " << pimod / lattice.volume()
                      << " ";
    }

    hila::out0 << "Moduli done \n";
}

void scaling_sim::write_energies() {
  
  Complex<double> suma(0),sumb2(0),sumb3(0),sumb4(0),sumb5(0);
  Complex<double> suma_we(0),sumb2_we(0),sumb3_we(0),sumb4_we(0),sumb5_we(0);
  Complex<double> sumk1(0),sumk2(0),sumk3(0);
  Complex<double> sumk1_we(0),sumk2_we(0),sumk3_we(0);
  real_t sumb1 = 0;
  Complex<double> sumb1_we = 0;
  Complex<double> sumkin(0);
  Complex<double> sumkin_we(0);


    
    hila::set_allreduce(false);
    onsites(ALL) {

      Matep MP;
      Complex<double> a(0),b2(0),b3(0),b4(0),b5(0);
      Complex<double> kin(0);
      Complex<double> k1(0), k2(0), k3(0);
      Complex<double> bfe(0);
      double b1 = 0;

      real_t ebfe=fmin(MP.f_A_td(p[X], T[X]),MP.f_B_td(p[X], T[X])); 
      
      real_t beta[6];
      point_params(T[X], p[X],beta);
      
      a = beta[0] * (A[X]*A[X].dagger()).trace();

      b1 = beta[1] * ((A[X]*A[X].transpose()).trace()).squarenorm();

      b2 = beta[2] * ((A[X]*A[X].dagger()).trace()*(A[X]*A[X].dagger()).trace());

      b3 = beta[3] * ((A[X]*A[X].transpose()*A[X].conj()*A[X].dagger()).trace());

      b4 = beta[4] * ((A[X]*A[X].dagger()*A[X]*A[X].dagger()).trace());

      b5 = beta[5] * ((A[X]*A[X].dagger()*A[X].conj()*A[X].transpose()).trace());

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
	              << sumb5_we.re / vol << " " << sumb5_we.im / vol << " ";
    }

       hila::out0 << "Energy done \n";

}

void scaling_sim::write_phases() {

  real_t ph0(0),ph1(0),ph2(0),ph3(0),ph4(0),ph5(0),ph6(0),ph7(0),ph8(0),ph9(0);

  hila::set_allreduce(false);

  onsites (ALL) {

    real_t R1,R2,R3,R4,R5;
    real_t p1,p2,p3,p4,p5,p6,p7,p8;
    real_t error=1.0/5.0;
    int phase;
    real_t p;
    phi_t Ac;

    Ac=A[X]/sqrt((A[X]*A[X].dagger()).trace());;

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

  hila::out0 << "Phases done \n";
  
}
void scaling_sim::write_positions() {

    
  hila::set_allreduce(false);

  if (config.write_phases==1)
    {
      Field<Vector<8,float>> data;
      
      onsites (ALL) {

	real_t R1,R2,R3,R4,R5;
	phi_t Ac;

	Ac=A[X]/sqrt((A[X]*A[X].dagger()).trace());;

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

	(A[X].dagger()*A[X]).eigen_jacobi(eval[X],evec[X]);

      }

      eval.write("points/eigenvalues-t"+std::to_string(int(t/config.dt)),false);
      evec.write("points/eigenvectors-t"+std::to_string(int(t/config.dt)),false);
      
    }
      
}

void scaling_sim::write_A_matrix_positions() {

   std::ofstream stream_out;
  
  const std::string fname = "A_matrix_output/t"+std::to_string(int(t/config.dt))+".dat";
  stream_out.open(fname, std::ios::out);

  hila::set_allreduce(false);
  A.write(stream_out,false,8);
  
  stream_out.close();

  // ******************************************************************************** //

}

void scaling_sim::latticeCoordinate_output()
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

void scaling_sim::next() {
   
  static hila::timer next_timer("timestep");
  Field<phi_t> deltaPi;
  Field<Vector<3,Complex<real_t>>> djAaj;

  int bc=config.boundary_conditions;
  
  next_timer.start();

  onsites (ALL) {

    Matep MP;
    real_t gapa = MP.gap_A_td(p[X], T[X]);
    real_t gapb = MP.gap_B_td(p[X], T[X]);

    
    A[X] += config.dt * pi[X];

    if (bc ==1)
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
    else if (bc==2)
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

    real_t beta[6];
    point_params(T[X], p[X],beta);
      
    auto AxAt = A[X]*A[X].transpose();
    auto AxAd = A[X]*A[X].dagger();
    
    deltaPi[X] = - beta[0]*A[X]
      - 2.0*beta[1]*A[X].conj()*AxAt.trace() 
      - 2.0*beta[2]*A[X]*AxAd.trace()
      - 2.0*beta[3]*AxAt*A[X].conj()
      - 2.0*beta[4]*AxAd*A[X] 
      - 2.0*beta[5]*A[X].conj()*A[X].transpose()*A[X];

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

    deltaPi[X] +=  (1.0 / (config.dx * config.dx)) * (A[X + e_x] + A[X - e_x]
						      + A[X + e_y] + A[X - e_y]
						      + A[X + e_z] + A[X - e_z]
						      - 6.0 * A[X]);
  }


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

void scaling_sim::next_bath() {

  static hila::timer next_timer("timestep");
  Field<phi_t> deltaPi;
  Field<Vector<3,Complex<real_t>>> djAaj;

  real_t ep2 = 1.0-exp(-2.0*config.gamma*config.dt) ;

  int bc=config.boundary_conditions;

    next_timer.start();

  onsites (ALL) {

    Matep MP;
    real_t gapa = MP.gap_A_td(p[X], T[X]);
    real_t gapb = MP.gap_B_td(p[X], T[X]);

    A[X] += config.dt * pi[X];

    if (bc ==1)
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
    else if (bc==2)
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

    real_t beta[6];
    point_params(T[X], p[X],beta);

    auto AxAt = A[X]*A[X].transpose();
    auto AxAd = A[X]*A[X].dagger();

    deltaPi[X] = - beta[0]*A[X]
      - 2.0*beta[1]*A[X].conj()*AxAt.trace()
      - 2.0*beta[2]*A[X]*AxAd.trace()
      - 2.0*beta[3]*AxAt*A[X].conj()
      - 2.0*beta[4]*AxAd*A[X]
      - 2.0*beta[5]*A[X].conj()*A[X].transpose()*A[X];

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

    real_t tb = T[X]/tc;
    real_t sig = sqrt(2.0*tb*config.gamma);
      
    deltaPi[X] += (1.0/(4.0*config.dx*config.dx)) * (A[X + e_x] + A[X - e_x]
                                                     + A[X + e_y] + A[X - e_y]
                                                     + A[X + e_z] + A[X - e_z]
                                                     - 6.0*A[X]);
     deltaPi[X] += hila::gaussrand()*sig;

  }

    //onsites (ALL) {deltaPi[X] *= config.dt;} // I think that this is the problem, multiplication with respect to dt                                                                                      \
                                                                                                                                                                                                            

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

	real_t beta[6];
	point_params(T[X], p[X],beta);
	
        auto AxAt = A[X]*A[X].transpose();
        auto AxAd = A[X]*A[X].dagger();

        deltaPi[X] = - beta[0]*A[X]
          - 2.0*beta[1]*A[X].conj()*AxAt.trace()
          - 2.0*beta[2]*A[X]*AxAd.trace()
          - 2.0*beta[3]*AxAt*A[X].conj()
          - 2.0*beta[4]*AxAd*A[X]
          - 2.0*beta[5]*A[X].conj()*A[X].transpose()*A[X];
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

	real_t tb = T[X]/tc;
	real_t sig = sqrt(2.0*tb*config.gamma);
	
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

void scaling_sim::nextT() {

  real_t te = t/config.tStart;
  
  onsites(ALL) {

    auto xcoord = X.coordinate(e_x);
    auto ycoord = X.coordinate(e_y);
    auto zcoord = X.coordinate(e_z);

    real_t r = sqrt((xcoord-config.lx/2.0)*(xcoord-config.lx/2.0)+(ycoord-config.ly/2.0)*(ycoord-config.ly/2.0)+(zcoord-config.lz/2.0)*(zcoord-config.lz/2.0));

    real_t expr=exp(-1.0*r*r/(config.diffT*te));

    T[X] = config.IniT + config.ampT*expr*pow(te,-3.0/2.0); //this is okay if diffT=sigT*sigT
    
  }
  
}


int main(int argc, char **argv) {
    scaling_sim sim;
    const std::string output_fname = sim.allocate("sim_params.txt", argc, argv);
    sim.initialize();
    sim.initializeT();
    sim.initializep();

    int stepspos;
    
    int steps =
        (sim.config.tEnd - sim.config.tStats) /
        (sim.config.dt * sim.config.nOutputs); // number of steps between printing stats
    if (steps == 0)
        steps = 1;

    if (sim.config.positions == 1) {
      stepspos = (sim.config.tEnd - sim.config.tStats) /
        (sim.config.dt * sim.config.npositionout);
      if (stepspos == 0)
        stepspos = 1;
    }
        
    int stat_counter = 0;

    if (hila::myrank() == 0) {
        sim.config.stream.open(output_fname, std::ios::out);
    }

      
    // on gpu the simulation timer is fake, because there's no sync here.  
    // BUt we want to avoid unnecessary sync anyway.
    static hila::timer run_timer("Simulation time"), meas_timer("Measurements");
    run_timer.start();

    // lattice coordinates output
    sim.latticeCoordinate_output();

    // output steps just before while loop
    // std::cout << "\n" << "steps is " << steps << " with nOutput " << sim.config.nOutputs << "\n";
    
    //auto tildephi = sim.phi;
    while (sim.t < sim.config.tEnd) {
        if (sim.t >= sim.config.tStats) {

	    // hila::out0 << "Writing output at time " << sim.t
	    //            << "stat_counter is " << stat_counter
            //            << "t point is " << sim.t/sim.config.dt
	    // 	       << "steps is " << steps
	    // 	       << "\n";
            
	  
            if (stat_counter % steps == 0) {
	      meas_timer.start();
	      hila::out0 << "Writing output at time " << sim.t 
			 << ", stat_counter is " << stat_counter
                         << ", t point is " << sim.t/sim.config.dt
			 << ", steps is " << steps
			 << "\n";
	      sim.write_moduli();
	      sim.write_energies();
	      sim.write_phases();
	      meas_timer.stop();
            }
	    if (sim.config.positions == 1)
	      {
		if (stat_counter % stepspos == 0)
		  {
		    sim.write_positions();
		  }
	      }
            stat_counter++;
        }
	if (sim.config.useTbath == 1)
	  {
	    sim.next_bath();}
	else
	  {
	    sim.next();}
	if(sim.config.evolveT == 1)
	  {
	    sim.nextT();
	  }
    }
    run_timer.stop();

    if (hila::myrank() == 0) {
        sim.config.stream.close();
    }

    hila::finishrun();
    return 0;
}

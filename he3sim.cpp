#define _USE_MATH_DEFINES
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
//#include <math.h>
#include <assert.h>

#include "plumbing/hila.h"
#include "plumbing/fft.h"


// Definition of the fiedl that we will use
using real_t = float;   // or double
using phi_t = Matrix<3,3,Complex<real_t>>;     // saves the trouble of writing this every time



/// Container for simulation parameters and methods
class scaling_sim{

public:
  scaling_sim() = default;
  const std::string allocate(const std::string &fname, int argc, char **argv);
  void initialize();
  void write_moduli();
  void write_energies();
  void next();
  
  Field<phi_t> A;
  Field<phi_t> pi;

  real_t t;
  
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
      
      real_t tau;
      real_t alpha;
      real_t beta1;
      real_t beta2;
      real_t beta3;
      real_t beta4;
      real_t beta5;

      real_t tStats;
      real_t nOutputs;
      std::fstream stream;
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
    config.tau = parameters.get("tau");
    config.alpha = parameters.get("alpha");
    config.beta1 = parameters.get("beta1");
    config.beta2 = parameters.get("beta2");
    config.beta3 = parameters.get("beta3");
    config.beta4 = parameters.get("beta4");
    config.beta5 = parameters.get("beta5");
    config.tStats = parameters.get("tStats");
    config.nOutputs = parameters.get("nOutputs");
    const std::string output_file = parameters.get("output_file");

    config.dt = config.dx * config.dtdxRatio;
    t = config.tStart;

    CoordinateVector box_dimensions = {config.l, config.l, config.l};
    lattice->setup(box_dimensions);
    hila::seed_random(config.seed);

    return output_file;
}


void scaling_sim::initialize() {

  int N = config.l;
  real_t dx = config.dx;

  output0<<"R: "<<hila::random()<<"\n";
  
  switch (config.initialCondition) {
    
  case 1: {
    pi = 0;
    foralldir(d1)foralldir(d2){
      onsites (ALL) {

	//if (d1==d2){ A[X].e(d1,d2)=1.0;}
	if(d1==d2){A[X].e(d1,d2).re = hila::random();}
	else {A[X].e(d1,d2).re = 0.0;}
	A[X].e(d1,d2).im = 0.0;//hila::random();

      }}

    onsites (ALL) {
      A[X] /= A[X].norm();
      A[X] *= 5.4772; //apparently the fix point
    }
    
        output0 << "Components randomly created \n";

        break;
    }
    
  default: {

    // #pragma hila ast_dump
    pi = 0.0; //set derivative matrix to zero
    onsites (ALL) {
      A[X].fill(1.0);
    }
    
    output0 << "Field matrix set to 1 everywhere \n";

    break;
  }
  }
}

void scaling_sim::write_moduli() {

   // real_t a = scaleFactor(t);

    double Amod = 0.0;
    double pimod = 0.0;

    hila::set_allreduce(false);
    onsites (ALL) {
        Amod += A[X].norm();
        pimod += pi[X].norm();
    }

    if (hila::myrank() == 0) {
        config.stream << t << " "
                      << Amod / lattice->volume() << " " << pimod / lattice->volume()
                      << " ";
    }
}

void scaling_sim::write_energies() {

  double sumar = 0.0;
  double sumai = 0.0;
  double sumb1 = 0.0;
  double sumb2r = 0.0;
  double sumb2i = 0.0;
  double sumb3r = 0.0;
  double sumb3i = 0.0;
  double sumb4r = 0.0;
  double sumb4i = 0.0;
  double sumb5r = 0.0;
  double sumb5i = 0.0;

  hila::set_allreduce(false);
  onsites (ALL) {

    sumar += config.alpha * ((A[X]*A[X].dagger()).trace()).re;
    sumai += config.alpha * ((A[X]*A[X].dagger()).trace()).im;
    
    sumb1 += config.beta1 * ((A[X]*A[X].transpose()).trace()).squarenorm();

    sumb2r += config.beta2 * ((A[X]*A[X].dagger()).trace()*(A[X]*A[X].dagger()).trace()).re;
    sumb2i += config.beta2 * ((A[X]*A[X].dagger()).trace()*(A[X]*A[X].dagger()).trace()).im;

    sumb3r += config.beta3 * ((A[X]*A[X].transpose()*A[X].conj()*A[X].dagger()).trace()).re;
    sumb3i += config.beta3 * ((A[X]*A[X].transpose()*A[X].conj()*A[X].dagger()).trace()).im;

    sumb4r += config.beta4 * ((A[X]*A[X].dagger()*A[X]*A[X].dagger()).trace()).re;
    sumb4i += config.beta4 * ((A[X]*A[X].dagger()*A[X]*A[X].dagger()).trace()).im;

    sumb5r += config.beta5 * ((A[X]*A[X].dagger()*A[X].conj()*A[X].transpose()).trace()).re;
    sumb5i += config.beta5 * ((A[X]*A[X].dagger()*A[X].conj()*A[X].transpose()).trace()).im;
  }

  if (hila::myrank() == 0) {
    double vol = (double)config.l * config.l * config.l;
    config.stream << sumar / vol << " "<< sumai / vol << " "
		  << sumb1 / vol << " "
		  << sumb2r / vol << " " << sumb2i / vol << " "
		  << sumb3r / vol << " " << sumb3i / vol << " "
 		  << sumb4r / vol << " " << sumb4i / vol << " "
		  << sumb5r / vol << " " << sumb5i / vol << "\n";
  }
}  
void scaling_sim::next() {


    static hila::timer next_timer("timestep");
    Field<phi_t> deltaPi;
    Field<Vector<3,Complex<real_t>>> djAaj;

    
    next_timer.start();

    onsites (ALL) {
        A[X] += config.dt * pi[X];
        deltaPi[X] = - config.alpha*A[X] - 2.0*config.beta1*A[X]*(A[X]*A[X].transpose()).trace() - 2.0*config.beta2*A[X]*(A[X]*A[X].dagger()).trace()
	  - 2.0*config.beta3*(A[X]*A[X].transpose()*A[X]) - 2.0*config.beta4*(A[X]*A[X].dagger()*A[X]) - 2.0*config.beta5*(A[X].conj()*A[X].transpose()*A[X]);
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
        deltaPi[X] += (1.0/(2.0*config.dx*config.dx))*mat;
    }

    onsites (ALL) {  
      deltaPi[X] += (1.0/(config.dx*config.dx)) * (A[X + e_x] + A[X -e_x]
						   + A[X + e_y] + A[X -e_y]
						   + A[X + e_z] + A[X -e_z]
						   - 6.0*A[X]);
    }
    

    if (t < config.tdif)
    {
        pi[ALL] = deltaPi[X]/(config.difFac*config.dt);
        t += config.dt/config.difFac;
    }
    else if (t < config.tdis && config.gamma > 0 )
    {
      pi[ALL] = pi[X] + (1.0/config.tau) * (deltaPi[X] - 2.0 * config.gamma * pi[X]);
      t += config.dt;
    }
    else
    {
      pi[ALL] = pi[X] + (1.0/config.tau) * deltaPi[X];
      t += config.dt;
    }

    next_timer.stop();

}

int main(int argc, char **argv) {
    scaling_sim sim;
    const std::string output_fname = sim.allocate("sim_params.txt", argc, argv);
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


    // on gpu the simulation timer is fake, because there's no sync here.  
    // BUt we want to avoid unnecessary sync anyway.
    static hila::timer run_timer("Simulation time"), meas_timer("Measurements");
    run_timer.start();
    
    //auto tildephi = sim.phi;
    while (sim.t < sim.config.tEnd) {
        if (sim.t >= sim.config.tStats) {
            if (stat_counter % steps == 0) {
	      meas_timer.start();
	      output0 << "Writing output at time " << sim.t << "\n";
	      sim.write_moduli();
	      sim.write_energies();
	      meas_timer.stop();
            }
            stat_counter++;
        }
        sim.next();
    }
    run_timer.stop();

    if (hila::myrank() == 0) {
        sim.config.stream.close();
    }

    hila::finishrun();
    return 0;
}

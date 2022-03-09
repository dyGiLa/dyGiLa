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

  switch (config.initialCondition) {
    
  case 1: {
    pi = 0;
    foralldir(d1) foralldir(d2){
      onsites (ALL) {

	A[X].e(d1,d2).re = hila::random();
	A[X].e(d1,d2).im = hila::random();
	
      }}

    onsites (ALL) {
      A[X] /= A[X].norm();}
    
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

  double suma = 0.0;
  double sumb1 = 0.0;
  double sumb2 = 0.0;
  double sumb3 = 0.0;
  double sumb4 = 0.0;
  double sumb5 = 0.0;

  hila::set_allreduce(false);
  onsites (ALL) {

    suma += config.alpha * ((A[X]*A[X].dagger()).trace()).squarenorm();
    sumb1 += config.beta1 * ((A[X]*A[X].transpose()).trace()).squarenorm();
    sumb2 += config.beta2 * ((A[X]*A[X].dagger()).trace()*(A[X]*A[X].dagger()).trace()).squarenorm();
    sumb3 += config.beta3 * ((A[X]*A[X].transpose()*A[X].conj()*A[X].dagger()).trace()).squarenorm();
    sumb4 += config.beta4 * ((A[X]*A[X].dagger()*A[X]*A[X].dagger()).trace()).squarenorm();
    sumb5 += config.beta5 * ((A[X]*A[X].dagger()*A[X].conj()*A[X].transpose()).trace()).squarenorm();
  }

  if (hila::myrank() == 0) {
    double vol = (double)config.l * config.l * config.l;
    config.stream << suma / vol << " " << sumb1 / vol << " " << sumb2 / vol << " " << sumb3 / vol << " " << sumb4 / vol << " " << sumb5 / vol << "\n ";
  }
}  
void scaling_sim::next() {


    static hila::timer next_timer("timestep");
    Field<phi_t> deltaPi;
    Field<phi_t> partialA;
    Field<phi_t> Appx;
    Field<phi_t> Appy;
    Field<phi_t> Appz;
    Field<phi_t> Ampx;
    Field<phi_t> Ampy;
    Field<phi_t> Ampz;
    Field<phi_t> Apmx;
    Field<phi_t> Apmy;
    Field<phi_t> Apmz;
    Field<phi_t> Ammx;
    Field<phi_t> Ammy;
    Field<phi_t> Ammz;
    
    next_timer.start();

    onsites (ALL) {
        A[X] += config.dt * pi[X];
        deltaPi[X] = - config.alpha*A[X] - 2.0*config.beta1*A[X]*(A[X]*A[X].transpose()).trace() + 2.0*config.beta2*A[X]*(A[X]*A[X].dagger()).trace()
	  + 2.0*config.beta3*(A[X]*A[X].transpose()*A[X]) + 2.0*config.beta4*(A[X]*A[X].dagger()*A[X]) + 2.0*config.beta5*(A[X].conj()*A[X].transpose()*A[X]);
    }

    foralldir(d1) foralldir(d2){
    onsites (ALL) {
      Appx[X]=A[X+d2+e_x];
      Appy[X]=A[X+d2+e_y];
      Appz[X]=A[X+d2+e_z];
      Ampx[X]=A[X-d2+e_x];
      Ampy[X]=A[X-d2+e_y];
      Ampz[X]=A[X-d2+e_z];
      Apmx[X]=A[X+d2-e_x];
      Apmy[X]=A[X+d2-e_y];
      Apmz[X]=A[X+d2-e_z];
      Ammx[X]=A[X-d2-e_x];
      Ammy[X]=A[X-d2-e_y];
      Ammz[X]=A[X-d2-e_z];
      
      deltaPi[X].e(d1,d2) += (1.0/(4.0*config.dx*config.dx)) * (Appx[X].e(d1,e_x) - Apmx[X].e(d1,e_x)
								+Appy[X].e(d1,e_y) - Apmy[X].e(d1,e_y)
								+Appz[X].e(d1,e_z) - Apmz[X].e(d1,e_z)
								-Ampx[X].e(d1,e_x) + Ammx[X].e(d1,e_x)
								-Ampy[X].e(d1,e_y) + Ammy[X].e(d1,e_y)
							        -Ampz[X].e(d1,e_z) + Ammz[X].e(d1,e_z));
    }}

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

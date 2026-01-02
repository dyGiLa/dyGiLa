#define USE_MPI 
// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <string>
// #include <assert.h>

#include "plumbing/hila.h"

#include "glsol.hpp"
#include "matep.hpp"


void glsol::dampAndRelax() {

  if (
      t < config.tdis
      && config.useTbath == 1
      && (extinguish_t * config.dt) < config.extinguish_off_t_count
     )
    {
      // though config.useTbath == 1 is till be true, only a big damping is needed.
      pi[ALL] = pi[X] + (deltaPi[X] - 1.0 * config.gamma * pi[X]) * config.dt;
      t += config.dt;
    }
  else if (
	   t < config.tdis
	   && config.useTbath == 1
	   && (extinguish_t * config.dt) >= config.extinguish_off_t_count
	  )
    {
      // though config.useTbath == 1 is till be true, no thermal noise is added, we just removed them.
      onsites(ALL)
	{
         matep::Matep MPonsites;
         pi[X] = pi[X] + (deltaPi[X] - 1.0 * MPonsites.gamma_td(config.Inip, T[X], phaseMarker[X]) * pi[X]) * config.dt;
	}

      t += config.dt;
    }  
  else
    {
      pi[ALL] = pi[X] + deltaPi[X]*config.dt;
      t += config.dt;
    }

} // dampAndRelax() ends here


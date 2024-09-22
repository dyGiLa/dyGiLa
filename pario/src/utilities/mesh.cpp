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

#include "ascent.hpp"
#include "conduit_blueprint.hpp"


void parIO::describeMesh(glsol &sol) {

    // Create a 3D mesh defined on a uniform grid of points
   
    // conduit::Node mesh;
    mesh["state/time"].set_external(&sol.t);
    // mesh["state/cycle"].set_external(&step);
#if defined USE_MPI
    mesh["state/domain_id"] = lattice.mynode.rank;
#endif
    mesh["state/software"] = "dyGiLa";
    mesh["state/title"] = "TDGL-Langiven equations simulator";
    mesh["state/info"] = "Parallel IO stream and insitu rendering of data from dyGiLa";

    // create the coordinate set
    mesh["coordsets/coords/type"] = "uniform";
    mesh["coordsets/coords/dims/i"] = lattice.mynode.size[0] + 2;
    mesh["coordsets/coords/dims/j"] = lattice.mynode.size[1] + 2;
    mesh["coordsets/coords/dims/k"] = lattice.mynode.size[2] + 2;

    // add origin and spacing to the coordset (optional)
    mesh["coordsets/coords/origin/x"] = ((lattice.mynode.min[0] - 1) * sol.config.dx);
    mesh["coordsets/coords/origin/y"] = ((lattice.mynode.min[1] - 1) * sol.config.dx);
    mesh["coordsets/coords/origin/z"] = ((lattice.mynode.min[2] - 1) * sol.config.dx);
    mesh["coordsets/coords/spacing/dx"] = sol.config.dx;
    mesh["coordsets/coords/spacing/dy"] = sol.config.dx;
    mesh["coordsets/coords/spacing/dz"] = sol.config.dx;
    
    // add the topology
    // this case is simple b/c it's implicitly derived from the coordinate set
    mesh["topologies/topo/type"] = "uniform";
    // reference the coordinate set by name
    mesh["topologies/topo/coordset"] = "coords";

} // describeMesh() end here


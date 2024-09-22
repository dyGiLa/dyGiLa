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


void parIO::describeMesh_addGhost_verify() {
    
    //create an element associated field named ghostCells
    mesh["fields/ascent_ghosts/association"] = "element";
    mesh["fields/ascent_ghosts/topology"] = "topo";
    mesh["fields/ascent_ghosts/values"].set_external(ghostCellsMask, ghostVolume);

    // make sure the mesh we created conforms to the blueprint
    conduit::Node verify_info;
    if (!conduit::blueprint::mesh::verify(mesh, verify_info)) {
        hila::out0 << "Mesh Verify failed!\n";
        hila::out0 << verify_info.to_yaml()
		   << '\n';
    }
    else {
        hila::out0 << "Mesh verify success!\n";
    }

} // describeMesh() end here


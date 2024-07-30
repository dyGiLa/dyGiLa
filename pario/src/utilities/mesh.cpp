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
#include "pario.hpp"

//#if defined USE_ASCENT

#include "ascent.hpp"
#include "conduit_blueprint.hpp"

//#endif

void parIO::describeMesh(glsol &sol) {

    // Create a 3D mesh defined on a uniform grid of points
   
    // conduit::Node mesh;
    mesh["state/time"].set_external(&sol.t);
    // mesh["state/cycle"].set_external(&step);
#if defined USE_MPI
    mesh["state/domain_id"] = lattice.mynode.rank;
#endif
    mesh["state/software"] = "dyGiLa";
    mesh["state/title"] = "TDGL equation simulator";
    mesh["state/info"] = "In Situ rendering of order parameter field from dyGiLa";

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

    // create an vertex associated field named gapAOrdered
    mesh["fields/gapAOrdered/association"] = "vertex";
    mesh["fields/gapAOrdered/topology"] = "topo";
    mesh["fields/gapAOrdered/values"].set_external(gapAOrdered.data(), latticeVolumeWithGhost);

    // create an vertex associated field named feDensityOrdered
    // mesh["fields/feDensityOrdered/association"] = "vertex";
    // mesh["fields/feDensityOrdered/topology"] = "topo";
    // mesh["fields/feDensityOrdered/values"].set_external(feDensityOrdered.data(), latticeVolumeWithGhost);

    // // create an vertex associated field named trAOrdered
    // if (config.hdf5_trA_output == 1){
    //   mesh["fields/trA_reOrdered/association"] = "vertex";
    //   mesh["fields/trA_reOrdered/topology"] = "topo";
    //   mesh["fields/trA_reOrdered/values"].set_external(trA_reOrdered.data(), latticeVolumeWithGhost);

    //   mesh["fields/trA_imOrdered/association"] = "vertex";
    //   mesh["fields/trA_imOrdered/topology"] = "topo";
    //   mesh["fields/trA_imOrdered/values"].set_external(trA_imOrdered.data(), latticeVolumeWithGhost);      
    // }

    // if (config.hdf5_eigvA_output == 1){
    //   mesh["fields/eigAv1Ordered/association"] = "vertex";
    //   mesh["fields/eigAv1Ordered/topology"] = "topo";
    //   mesh["fields/eigAv1Ordered/values"].set_external(eigAv1Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/eigAv2Ordered/association"] = "vertex";
    //   mesh["fields/eigAv2Ordered/topology"] = "topo";
    //   mesh["fields/eigAv2Ordered/values"].set_external(eigAv2Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/eigAv3Ordered/association"] = "vertex";
    //   mesh["fields/eigAv3Ordered/topology"] = "topo";
    //   mesh["fields/eigAv3Ordered/values"].set_external(eigAv3Ordered.data(), latticeVolumeWithGhost);      
    // }

    // if (config.hdf5_mass_current_output == 1){
    //   mesh["fields/jm1Ordered/association"] = "vertex";
    //   mesh["fields/jm1Ordered/topology"] = "topo";
    //   mesh["fields/jm1Ordered/values"].set_external(jm1Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/jm2Ordered/association"] = "vertex";
    //   mesh["fields/jm2Ordered/topology"] = "topo";
    //   mesh["fields/jm2Ordered/values"].set_external(jm2Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/jm3Ordered/association"] = "vertex";
    //   mesh["fields/jm3Ordered/topology"] = "topo";
    //   mesh["fields/jm3Ordered/values"].set_external(jm3Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/phaseExpModulusO/association"] = "vertex";
    //   mesh["fields/phaseExpModulusO/topology"] = "topo";
    //   mesh["fields/phaseExpModulusO/values"].set_external(phaseExpModulusO.data(), latticeVolumeWithGhost);

    //   mesh["fields/phaseExpAngleO/association"] = "vertex";
    //   mesh["fields/phaseExpAngleO/topology"] = "topo";
    //   mesh["fields/phaseExpAngleO/values"].set_external(phaseExpAngleO.data(), latticeVolumeWithGhost);

    //   mesh["fields/phaseExp2ReO/association"] = "vertex";
    //   mesh["fields/phaseExp2ReO/topology"] = "topo";
    //   mesh["fields/phaseExp2ReO/values"].set_external(phaseExp2ReO.data(), latticeVolumeWithGhost);

    //   mesh["fields/phaseExp2ImO/association"] = "vertex";
    //   mesh["fields/phaseExp2ImO/topology"] = "topo";
    //   mesh["fields/phaseExp2ImO/values"].set_external(phaseExp2ImO.data(), latticeVolumeWithGhost);      
    // }    

    // if (config.hdf5_spin_current_output == 1){
    //   mesh["fields/js11O/association"] = "vertex";
    //   mesh["fields/js11O/topology"] = "topo";
    //   mesh["fields/js11O/values"].set_external(js11O.data(), latticeVolumeWithGhost);

    //   mesh["fields/js21O/association"] = "vertex";
    //   mesh["fields/js21O/topology"] = "topo";
    //   mesh["fields/js21O/values"].set_external(js21O.data(), latticeVolumeWithGhost);

    //   mesh["fields/js31O/association"] = "vertex";
    //   mesh["fields/js31O/topology"] = "topo";
    //   mesh["fields/js31O/values"].set_external(js31O.data(), latticeVolumeWithGhost);

    //   mesh["fields/js12O/association"] = "vertex";
    //   mesh["fields/js12O/topology"] = "topo";
    //   mesh["fields/js12O/values"].set_external(js12O.data(), latticeVolumeWithGhost);

    //   mesh["fields/js22O/association"] = "vertex";
    //   mesh["fields/js22O/topology"] = "topo";
    //   mesh["fields/js22O/values"].set_external(js22O.data(), latticeVolumeWithGhost);

    //   mesh["fields/js32O/association"] = "vertex";
    //   mesh["fields/js32O/topology"] = "topo";
    //   mesh["fields/js32O/values"].set_external(js32O.data(), latticeVolumeWithGhost);

    //   mesh["fields/js13O/association"] = "vertex";
    //   mesh["fields/js13O/topology"] = "topo";
    //   mesh["fields/js13O/values"].set_external(js13O.data(), latticeVolumeWithGhost);

    //   mesh["fields/js23O/association"] = "vertex";
    //   mesh["fields/js23O/topology"] = "topo";
    //   mesh["fields/js23O/values"].set_external(js23O.data(), latticeVolumeWithGhost);

    //   mesh["fields/js33O/association"] = "vertex";
    //   mesh["fields/js33O/topology"] = "topo";
    //   mesh["fields/js33O/values"].set_external(js33O.data(), latticeVolumeWithGhost);
   
    // }    

    
    // /*----------------------------------------------------------------------*/
    // /*---create vertices associated field named of uxxOrdered vxxOrdered ---*/
    // /*----------------------------------------------------------------------*/
    // if (config.A_matrix_output == 1){
    //   mesh["fields/u11Ordered/association"] = "vertex";
    //   mesh["fields/u11Ordered/topology"] = "topo";
    //   mesh["fields/u11Ordered/values"].set_external(u11Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/u12Ordered/association"] = "vertex";
    //   mesh["fields/u12Ordered/topology"] = "topo";
    //   mesh["fields/u12Ordered/values"].set_external(u12Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/u13Ordered/association"] = "vertex";
    //   mesh["fields/u13Ordered/topology"] = "topo";
    //   mesh["fields/u13Ordered/values"].set_external(u13Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/u21Ordered/association"] = "vertex";
    //   mesh["fields/u21Ordered/topology"] = "topo";
    //   mesh["fields/u21Ordered/values"].set_external(u21Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/u22Ordered/association"] = "vertex";
    //   mesh["fields/u22Ordered/topology"] = "topo";
    //   mesh["fields/u22Ordered/values"].set_external(u22Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/u23Ordered/association"] = "vertex";
    //   mesh["fields/u23Ordered/topology"] = "topo";
    //   mesh["fields/u23Ordered/values"].set_external(u23Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/u31Ordered/association"] = "vertex";
    //   mesh["fields/u31Ordered/topology"] = "topo";
    //   mesh["fields/u31Ordered/values"].set_external(u31Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/u32Ordered/association"] = "vertex";
    //   mesh["fields/u32Ordered/topology"] = "topo";
    //   mesh["fields/u32Ordered/values"].set_external(u32Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/u33Ordered/association"] = "vertex";
    //   mesh["fields/u33Ordered/topology"] = "topo";
    //   mesh["fields/u33Ordered/values"].set_external(u33Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/v11Ordered/association"] = "vertex";
    //   mesh["fields/v11Ordered/topology"] = "topo";
    //   mesh["fields/v11Ordered/values"].set_external(v11Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/v12Ordered/association"] = "vertex";
    //   mesh["fields/v12Ordered/topology"] = "topo";
    //   mesh["fields/v12Ordered/values"].set_external(v12Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/v13Ordered/association"] = "vertex";
    //   mesh["fields/v13Ordered/topology"] = "topo";
    //   mesh["fields/v13Ordered/values"].set_external(v13Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/v21Ordered/association"] = "vertex";
    //   mesh["fields/v21Ordered/topology"] = "topo";
    //   mesh["fields/v21Ordered/values"].set_external(v21Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/v22Ordered/association"] = "vertex";
    //   mesh["fields/v22Ordered/topology"] = "topo";
    //   mesh["fields/v22Ordered/values"].set_external(v22Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/v23Ordered/association"] = "vertex";
    //   mesh["fields/v23Ordered/topology"] = "topo";
    //   mesh["fields/v23Ordered/values"].set_external(v23Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/v31Ordered/association"] = "vertex";
    //   mesh["fields/v31Ordered/topology"] = "topo";
    //   mesh["fields/v31Ordered/values"].set_external(v31Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/v32Ordered/association"] = "vertex";
    //   mesh["fields/v32Ordered/topology"] = "topo";
    //   mesh["fields/v32Ordered/values"].set_external(v32Ordered.data(), latticeVolumeWithGhost);

    //   mesh["fields/v33Ordered/association"] = "vertex";
    //   mesh["fields/v33Ordered/topology"] = "topo";
    //   mesh["fields/v33Ordered/values"].set_external(v33Ordered.data(), latticeVolumeWithGhost);              
    // }
    /*----------------------------------------------------------------------*/
    /*--- vertex associated field named of uxxOrdered vxxOrdered end here --*/
    /*----------------------------------------------------------------------*/

    
    // create an element associated field named ghostCells
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


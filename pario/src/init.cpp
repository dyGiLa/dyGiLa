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

void parIO::init(glsol &sol) {

    latticeVolumeWithGhost =
        (lattice.mynode.size[0] + 2) * (lattice.mynode.size[1] + 2) * (lattice.mynode.size[2] + 2);
    
    latticeVolume =
        (lattice.mynode.size[0]) * (lattice.mynode.size[1]) * (lattice.mynode.size[2]);

    gapAOrdered.reserve(latticeVolumeWithGhost);
    feDensityOrdered.reserve(latticeVolumeWithGhost);
    Temperature.reserve(latticeVolumeWithGhost);
    phaseMarker.reserve(latticeVolumeWithGhost);

    // if (config.hdf5_trA_output == 1){
    //   trA_reOrdered.reserve(latticeVolumeWithGhost);
    //   trA_imOrdered.reserve(latticeVolumeWithGhost);      
    // }

    // if (config.hdf5_eigvA_output == 1){
    //   eigAv1Ordered.reserve(latticeVolumeWithGhost);
    //   eigAv2Ordered.reserve(latticeVolumeWithGhost);
    //   eigAv3Ordered.reserve(latticeVolumeWithGhost);      
    // }

    if (sol.config.hdf5_mass_current_output == 1){
      jm1Ordered.reserve(latticeVolumeWithGhost);
      jm2Ordered.reserve(latticeVolumeWithGhost);
      jm3Ordered.reserve(latticeVolumeWithGhost);

      phaseExpModulusO.reserve(latticeVolumeWithGhost);
      phaseExpAngleO.reserve(latticeVolumeWithGhost);
      phaseExp2ReO.reserve(latticeVolumeWithGhost);
      phaseExp2ImO.reserve(latticeVolumeWithGhost);                  
    }

    if (sol.config.hdf5_spin_current_output == 1){
      js11O.reserve(latticeVolumeWithGhost);
      js21O.reserve(latticeVolumeWithGhost);
      js31O.reserve(latticeVolumeWithGhost);

      js12O.reserve(latticeVolumeWithGhost);
      js22O.reserve(latticeVolumeWithGhost);
      js32O.reserve(latticeVolumeWithGhost);

      js13O.reserve(latticeVolumeWithGhost);
      js23O.reserve(latticeVolumeWithGhost);
      js33O.reserve(latticeVolumeWithGhost);      
    }

    if (sol.config.hdf5_A_matrix_output == 1){
     u11Ordered.reserve(latticeVolumeWithGhost); v11Ordered.reserve(latticeVolumeWithGhost);
     u12Ordered.reserve(latticeVolumeWithGhost); v12Ordered.reserve(latticeVolumeWithGhost);
     u13Ordered.reserve(latticeVolumeWithGhost); v13Ordered.reserve(latticeVolumeWithGhost);
     u21Ordered.reserve(latticeVolumeWithGhost); v21Ordered.reserve(latticeVolumeWithGhost);
     u22Ordered.reserve(latticeVolumeWithGhost); v22Ordered.reserve(latticeVolumeWithGhost);
     u23Ordered.reserve(latticeVolumeWithGhost); v23Ordered.reserve(latticeVolumeWithGhost);
     u31Ordered.reserve(latticeVolumeWithGhost); v31Ordered.reserve(latticeVolumeWithGhost);
     u32Ordered.reserve(latticeVolumeWithGhost); v32Ordered.reserve(latticeVolumeWithGhost);
     u33Ordered.reserve(latticeVolumeWithGhost); v33Ordered.reserve(latticeVolumeWithGhost);    
    }
     
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

    /*********************************/
    /*    all describeMesh calls     */
    /*********************************/    
    describeMesh(sol);
    describeMesh_gapA_FEDensity();
    describeMesh_Temperature();
    describeMesh_phaseMarker();

    if (sol.config.hdf5_mass_current_output == 1) {describeMesh_gapA_FEDensity();}
    if (sol.config.hdf5_spin_current_output == 1) {describeMesh_massCurrent();}
    if (sol.config.hdf5_A_matrix_output == 1) {describeMesh_AMatrix();}

    describeMesh_addGhost_verify();
    /*********************************/
    /*  describeMesh calls end here  */
    /*********************************/    
    

    pio_options["mpi_comm"] = MPI_Comm_c2f(lattice.mpi_comm_lat);
    pio_options["runtime/type"] = "ascent";
#if defined CUDA
    pio_options["runtime/vtkm/backend"] = "cuda";
    pio_options["cuda/init"] = "false";
#endif
    pio_options["timings"] = "false";
    
    pio.open(pio_options);
    pio.publish(mesh);


    /*********************************/
    /*   all defineActions calls     */
    /*********************************/        
    defineActions_insitu(sol);

    if (sol.config.hdf5_mass_current_output == 1) {defineActions_massCurrent();}
    if (sol.config.hdf5_spin_current_output == 1) {defineActions_spinCurrent();}
    if (sol.config.hdf5_A_matrix_output == 1) {defineActions_AMatrix();}

    defineActions_printTree();
    
    /*********************************/
    /*   all defineActions calls     */
    /*********************************/        

    
} // init() end here


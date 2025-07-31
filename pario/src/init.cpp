#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
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

    gapAContainer.reserve(latticeVolumeWithGhost);
    
    if (sol.config.pario_compute_feDensity == 1) { feDensityContainer.reserve(latticeVolumeWithGhost); }
    
    Temperature.reserve(latticeVolumeWithGhost);
    phaseMarker.reserve(latticeVolumeWithGhost);

    if (sol.config.hdf5_mass_current_output == 1){
      jm1Container.reserve(latticeVolumeWithGhost);
      jm2Container.reserve(latticeVolumeWithGhost);
      jm3Container.reserve(latticeVolumeWithGhost);

      phaseExpModulusContainer.reserve(latticeVolumeWithGhost);
      phaseExpAngleContainer.reserve(latticeVolumeWithGhost);
      phaseExp2ReContainer.reserve(latticeVolumeWithGhost);
      phaseExp2ImContainer.reserve(latticeVolumeWithGhost);                  
    }

    if (sol.config.hdf5_spin_current_output == 1){
      js11Container.reserve(latticeVolumeWithGhost);
      js21Container.reserve(latticeVolumeWithGhost);
      js31Container.reserve(latticeVolumeWithGhost);

      js12Container.reserve(latticeVolumeWithGhost);
      js22Container.reserve(latticeVolumeWithGhost);
      js32Container.reserve(latticeVolumeWithGhost);

      js13Container.reserve(latticeVolumeWithGhost);
      js23Container.reserve(latticeVolumeWithGhost);
      js33Container.reserve(latticeVolumeWithGhost);      
    }

    if (sol.config.hdf5_A_matrix_output == 1){
     u11Container.reserve(latticeVolumeWithGhost); v11Container.reserve(latticeVolumeWithGhost);
     u12Container.reserve(latticeVolumeWithGhost); v12Container.reserve(latticeVolumeWithGhost);
     u13Container.reserve(latticeVolumeWithGhost); v13Container.reserve(latticeVolumeWithGhost);
     u21Container.reserve(latticeVolumeWithGhost); v21Container.reserve(latticeVolumeWithGhost);
     u22Container.reserve(latticeVolumeWithGhost); v22Container.reserve(latticeVolumeWithGhost);
     u23Container.reserve(latticeVolumeWithGhost); v23Container.reserve(latticeVolumeWithGhost);
     u31Container.reserve(latticeVolumeWithGhost); v31Container.reserve(latticeVolumeWithGhost);
     u32Container.reserve(latticeVolumeWithGhost); v32Container.reserve(latticeVolumeWithGhost);
     u33Container.reserve(latticeVolumeWithGhost); v33Container.reserve(latticeVolumeWithGhost);    
    }
     
    // Containerne more point in each direction, but cell data (Npts - 1 cells)
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
    describeMesh_gapA_FEDensity(sol);
    describeMesh_Temperature();
    describeMesh_phaseMarker();

    if ((!!sol.config.hdf5Ststart == true) && (!!sol.config.hdf5Stend == true))
      {
       if (sol.config.hdf5_mass_current_output == 1) {describeMesh_massCurrent();}
       if (sol.config.hdf5_spin_current_output == 1) {describeMesh_spinCurrent();}
       if (sol.config.hdf5_A_matrix_output == 1) {describeMesh_AMatrix();}
      }

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

    if ((!!sol.config.hdf5Ststart == true) && (!!sol.config.hdf5Stend == true))
      {
       if (sol.config.hdf5_mass_current_output == 1) {defineActions_massCurrent(sol);}
       if (sol.config.hdf5_spin_current_output == 1) {defineActions_spinCurrent(sol);}
       if (sol.config.hdf5_A_matrix_output == 1) {defineActions_AMatrix(sol);}
       if (sol.config.hdf5_pMarker_output == 1) {defineActions_phaseMarker(sol);}
      }
    

    defineActions_printTree();
    
    /*********************************/
    /*   all defineActions calls     */
    /*********************************/        

    
} // init() end here


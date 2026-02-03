#define USE_BOTHSIDE_GHOSTS
#define USE_ADGRZ
// #define USE_MPI 
// #include <sstream>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <string>
// #include <assert.h>

#include "plumbing/hila.h"
//#include "plumbing/fft.h"
#include "plumbing/memalloc.h"

#include "glsol.hpp"
//#include "matep.hpp"
#include "pario.hpp"

// #include "ascent.hpp"
// #include "conduit_blueprint.hpp"

void parIO::ghostMask(glsol &sol) {
    
    // cell data (Npts - 1 cells) in each directions
    auto ghostNX = lattice->mynode.size[0] + 2 - 1;
    auto ghostNY = lattice->mynode.size[1] + 2 - 1;
    auto ghostNZ = lattice->mynode.size[2] + 2 - 1;

#ifdef USE_ADGRZ
    //auto mynodeExtentZ = lattice->mynode.size[2]; 
    auto mynodeMinZcoord = lattice->mynode.min[2];
    auto mynodeMaxZcoord = lattice->mynode.min[2] + (lattice->mynode.size[2]-1 ); // z-roof coordinat index: min+(size-1) 

    // hila::out << "mynodeMinZcoord is " << mynodeMinZcoord
    // 	      << ", mynodeMaxZcoord is " << mynodeMaxZcoord
    // 	      << std::endl;
#endif    
    
    ghostVolume = ghostNX * ghostNY * ghostNZ;
    ghostCellsMask = (unsigned char *)memalloc(ghostVolume * sizeof(unsigned char));

    long long counter = 0;
    unsigned char Mask = 0;

    // i,j,k = 0 left-most halos
    for (auto k = 0; k < ghostNZ; k++) {
        for (auto j = 0; j < ghostNY; j++) {
            for (auto i = 0; i < ghostNX; i++) {
#if !defined(USE_ADGRZ) && defined(USE_BOTHSIDE_GHOSTS)
	      bool kGhostFlag = (k == 0) || (k == (ghostNZ -1));
              bool jGhostFlag = (j == 0) || (j == (ghostNY -1));
              bool iGhostFlag = (i == 0) || (i == (ghostNX -1));	      
#elif defined(USE_ADGRZ) && defined(USE_BOTHSIDE_GHOSTS)
              bool kGhostFlag = ((k == 0) && !(mynodeMinZcoord == 0))
		                || ((k == (ghostNZ -1)) && !(mynodeMaxZcoord == (sol.config.lz-1)));
              bool jGhostFlag = (j == 0) || (j == (ghostNY -1));
              bool iGhostFlag = (i == 0) || (i == (ghostNX -1));	      	      
#elif !defined(USE_ADGRZ) && !defined(USE_BOTHSIDE_GHOSTS)
	      bool kGhostFlag = (k == 0);
              bool jGhostFlag = (j == 0);
              bool iGhostFlag = (i == 0);	      	      
#elif defined(USE_ADGRZ) && !defined(USE_BOTHSIDE_GHOSTS)
              bool kGhostFlag = ((k == 0) && !(mynodeMinZcoord == 0));
              bool jGhostFlag = (j == 0);
              bool iGhostFlag = (i == 0);	      	      	      
#endif		
              Mask = static_cast<unsigned char>(iGhostFlag || jGhostFlag || kGhostFlag);
              ghostCellsMask[counter++] = Mask;
            }
        }
    }
    
} // ghostMask() end here


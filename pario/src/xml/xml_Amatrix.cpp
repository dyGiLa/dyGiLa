#define USE_ASCENT 
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


void parIO::xml_Amatrix(glsol &sol){

  const std::string fname = "rank_xmls/" + sol.config.xmf2_Amatrix_fname + "_" + std::to_string(hila::myrank()) + ".xml";
  xml_out.open(fname, std::ios::out);
  
  const long dim_0 = lattice.mynode.size[0] + 2,
             dim_1 = lattice.mynode.size[1] + 2,
             dim_2 = lattice.mynode.size[2] + 2;

  unsigned int n;

  xml_out << "<Grid Name=\"dyGiLa-sim-Amatrix\" Type=\"Uniform\">\n"
          << "  <Topology name=\"topo\" TopologyType=\"3DRectMesh\" Dimensions=\""
          << dim_2 << " " << dim_1 << " " << dim_0 << "\"" << ">" << "\n"
          << "  </Topology>\n"
          << "  <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n"
          << "    <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">"
	  << "\n"
          << "      "
          << ((lattice.mynode.min[0] - 1) * sol.config.dx) << " "
          << ((lattice.mynode.min[1] - 1) * sol.config.dx) << " "
          << ((lattice.mynode.min[2] - 1) * sol.config.dx) << "\n"
          << "    </DataItem>\n"
          << "    <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n"
          << "      "
          << sol.config.dx << " " << sol.config.dx << " " << sol.config.dx << "\n"
          << "    </DataItem>\n"
          << "  </Geometry>\n"
          << "  <Attribute Name=\"gapA\" AttributeType=\"Scalar\" Center=\"Node\">\n"
          << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
          << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	  << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank() << ".hdf5:/mesh/fields/gapA/values"
	  << "\n"
          << "   </DataItem>\n"
          << "  </Attribute>\n"
          << "  <Attribute Name=\"feDensity\" AttributeType=\"Scalar\" Center=\"Node\">\n"
          << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
          << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	  << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank() << ".hdf5:/mesh/fields/feDensity/values"
	  << "\n"
          << "   </DataItem>\n"
          << "  </Attribute>"
 	  << "\n"
	  << std::flush;
    

   for (n = 0; n<=8; ++n){
         xml_out << "  <Attribute Name=\"u"
	         << std::to_string(n/3u + 1)
	      	 << std::to_string(n%3u + 1) << "\""
	         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                 << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                 << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
	         << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
	         << ".hdf5:/mesh/fields/u"
	         << std::to_string(n/3u + 1)
	         << std::to_string(n%3u + 1) << "Container/values"
	         << "\n"
                 << "   </DataItem>\n"
                 << "  </Attribute>"
	         << "\n"
	         << "  <Attribute Name=\"v"
	         << std::to_string(n/3u + 1)
	         << std::to_string(n%3u + 1) << "\""
	         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
                 << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                 << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
	         << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
	         << ".hdf5:/mesh/fields/v"
	         << std::to_string(n/3u + 1)
	         << std::to_string(n%3u + 1) << "Container/values"
	         << "\n"
                 << "   </DataItem>\n"
                 << "  </Attribute>"
	         << "\n" << std::flush;
   } // for n < 8 loop block
     
   xml_out << "  <Attribute Name=\"vtkGhostType\" AttributeType=\"Scalar\" Center=\"Cell\">\n"
           << "   <DataItem Format=\"HDF\" DataType=\"UChar\" Dimensions=\""
           << dim_0 - 1 << " " << dim_1 - 1 << " " << dim_2 - 1 << "\"" << ">" << "\n"
           << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank() << ".hdf5:/mesh/fields/ascent_ghosts/values"
           << "\n"   
           << "   </DataItem>\n"
           << "  </Attribute>\n"
           << "</Grid>"
           << "\n"
	   << std::flush;

  xml_out.close();

} // xml() end here


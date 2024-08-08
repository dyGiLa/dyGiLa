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

#include "ascent.hpp"
#include "conduit_blueprint.hpp"


void parIO::xml(glsol &sol){

  const std::string fname = "rank_xmls/" + sol.config.xmf2_fname + "_" + std::to_string(hila::myrank()) + ".xml";
  sol.config.xml_out.open(fname, std::ios::out);

  // hila::out << sol.config.xml_out.good() << "\n";
  // hila::out << sol.config.xml_out.is_open() << "\n";
  // hila::out << sol.config.xml_out.fail() << "\n";
  
  const long dim_0 = lattice.mynode.size[0] + 2,
             dim_1 = lattice.mynode.size[1] + 2,
             dim_2 = lattice.mynode.size[2] + 2;

  unsigned int n;

  sol.config.xml_out << "<Grid Name=\"sim-data\" Type=\"Uniform\">\n"
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
 	             << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank() << ".hdf5:/fields/gapAOrdered/values"
		     << "\n"
                     << "   </DataItem>\n"
                     << "  </Attribute>\n"
                     << "  <Attribute Name=\"feDensity\" AttributeType=\"Scalar\" Center=\"Node\">\n"
                     << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
                     << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
 	             << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank() << ".hdf5:/fields/feDensityOrdered/values"
		     << "\n"
                     << "   </DataItem>\n"
                     << "  </Attribute>"
 		     << "\n"
		     << std::flush;
    
    // if (config.hdf5_A_matrix_output == 1){
    //   for (n = 0; n<=8; ++n){
    //      config.xml_out << "  <Attribute Name=\"u"
    // 	                 << std::to_string(n/3u + 1)
    // 	         	 << std::to_string(n%3u + 1) << "\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/u"
    // 		         << std::to_string(n/3u + 1)
    // 		         << std::to_string(n%3u + 1) << "Ordered/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"
    // 	                 << "\n"
    // 	                 << "  <Attribute Name=\"v"
    // 	                 << std::to_string(n/3u + 1)
    // 		         << std::to_string(n%3u + 1) << "\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/v"
    // 		         << std::to_string(n/3u + 1)
    // 		         << std::to_string(n%3u + 1) << "Ordered/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"
    // 	                 << "\n" << std::flush;
    //   }

    // }

    // if (config.hdf5_trA_output == 1) {
    //      config.xml_out << "  <Attribute Name=\"trA_re\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/trA_reOrdered/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"
    // 	                 << "\n"
    // 	                 << "  <Attribute Name=\"trA_im\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/trA_imOrdered/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"
    // 	                 << "\n" << std::flush;
    // }

    // if (config.hdf5_eigvA_output == 1){
    //      config.xml_out << "  <Attribute Name=\"eigVal1_A\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/eigAv1Ordered/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"
    // 	                 << "\n"
    // 	                 << "  <Attribute Name=\"eigVal2_A\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/eigAv2Ordered/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"
    // 			 << "\n"  
    // 	                 << "  <Attribute Name=\"eigVal3_A\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/eigAv3Ordered/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"
    // 			 << "\n" << std::flush;
    // }

    // if (config.hdf5_mass_current_output == 1){
    //      config.xml_out << "  <Attribute Name=\"jm_1\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/jm1Ordered/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"
    // 	                 << "\n"
    // 	                 << "  <Attribute Name=\"jm_2\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/jm2Ordered/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"
    // 			 << "\n"  
    // 	                 << "  <Attribute Name=\"jm_3\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/jm3Ordered/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"
    // 	                 << "  <Attribute Name=\"phaseExpModulus\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/phaseExpModulusO/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"
    // 	                 << "  <Attribute Name=\"phaseExpAngle\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/phaseExpAngleO/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"
    // 	                 << "  <Attribute Name=\"phaseExp2ReO\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/phaseExp2ReO/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"
    // 	                 << "  <Attribute Name=\"phaseExp2ImO\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/phaseExp2ImO/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"	   	   	   
    // 			 << "\n" << std::flush;
    // }

    // if (config.hdf5_spin_current_output == 1){
    //       config.xml_out << "  <Attribute Name=\"js_11\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/js11O/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"
    // 	                 << "\n"
    // 	                 << "  <Attribute Name=\"js_21\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/js21O/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"
    // 			 << "\n"  
    // 	                 << "  <Attribute Name=\"js_31\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/js31O/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"
    //                      << "\n"	    
    //                      << "  <Attribute Name=\"js_12\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/js12O/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"
    // 	                 << "\n"
    // 	                 << "  <Attribute Name=\"js_22\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/js22O/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"
    // 			 << "\n"  
    // 	                 << "  <Attribute Name=\"js_32\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/js32O/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"	    
    //                      << "\n"
    //                      << "  <Attribute Name=\"js_13\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/js13O/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"
    // 	                 << "\n"
    // 	                 << "  <Attribute Name=\"js_23\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/js23O/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"
    // 			 << "\n"  
    // 	                 << "  <Attribute Name=\"js_33\""
    // 		         << " AttributeType=\"Scalar\" Center=\"Node\">\n"
    //                      << "   <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=" << "\""
    //                      << dim_0 << " " << dim_1 << " " << dim_2 << "\"" << ">" << "\n"
    // 	                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank()
    // 		         << ".hdf5:/fields/js33O/values"
    // 		         << "\n"
    //                      << "   </DataItem>\n"
    //                      << "  </Attribute>"	    
    // 			 << "\n" << std::flush;
    // }
     
  sol.config.xml_out << "  <Attribute Name=\"vtkGhostType\" AttributeType=\"Scalar\" Center=\"Cell\">\n"
                 << "   <DataItem Format=\"HDF\" DataType=\"UChar\" Dimensions=\""
                 << dim_0 - 1 << " " << dim_1 - 1 << " " << dim_2 - 1 << "\"" << ">" << "\n"
                 << "    domain_" << std::setfill('0') << std::setw(6) << hila::myrank() << ".hdf5:/fields/ascent_ghosts/values"
         	 << "\n"   
                 << "   </DataItem>\n"
                 << "  </Attribute>\n"
                 << "</Grid>"
                 << "\n"
		 << std::flush;

  sol.config.xml_out.close();

} // xml() end here


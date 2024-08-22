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

void parIO::xdmf(glsol &sol){

  unsigned int rank_no = 0/*hila::myrank()*/;
  std::fstream xml_file;
  xdmf_out.open(sol.config.xmf2_fname, std::ios::out);

  xdmf_out << "<?xml version=\"1.0\" ?>\n"
           << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
           << "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.0\">\n"
           << "\n"
           << "<Domain>\n"
           << "<Grid GridType=\"Collection\" CollectionType=\"Collection\">"
           << "\n";
  
  while(rank_no < hila::number_of_nodes()){
    const std::string fname = "rank_xmls/" + sol.config.xmf2_fname + "_" + std::to_string(rank_no) + ".xml";
    xml_file.open(fname, std::ios::in);
    xdmf_out << xml_file.rdbuf() << "\n" << std::endl;
    xml_file.close();

    ++rank_no;
  }

  xdmf_out << "</Grid>"
           << "\n"
           << "</Domain>\n"
           << "</Xdmf>"
           << std::endl;
  
  xdmf_out.close();
} // xdmf() end here


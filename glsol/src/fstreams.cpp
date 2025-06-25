#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>
#include <stdlib.h>

#include "glsol.hpp"

void glsol::fstreams_open(const std::vector<std::string> &name_files) {

  /* measure-stream file open */
  config.stream.open(name_files[0], std::ios::out);

  config.stream << "\"t\"" << "," << "\"T_000Q\"" << ","
                /***************************/
                << "\"sumAgapRe_VA\"" << "," << "\"sumAgapIm_VA\"" << ","
                << "\"sumgapARe_VA\"" << "," << "\"sumgapAIm_VA\"" << ","
                /***************************/
                << "\"sumkinRe_VA\"" << "," << "\"sumkinIm_VA\"" << ","
                << "\"sumkin_weRe_VA\"" << "," << "\"sumkin_weIm_VA\"" << ","    
                /***************************/
                << "\"sumk1Re_VA\"" << "," << "\"sumk1Im_VA\"" << ","
                << "\"sumk1_weRe_VA\"" << "," << "\"sumk1_weIm_VA\"" << ","    
                /***************************/
                << "\"sumk2Re_VA\"" << "," << "\"sumk2Im_VA\"" << ","
                << "\"sumk2_weRe_VA\"" << "," << "\"sumk2_weIm_VA\"" << ","    
                /***************************/
                << "\"sumk3Re_VA\"" << "," << "\"sumk3Im_VA\"" << ","
                << "\"sumk3_weRe_VA\"" << "," << "\"sumk3_weIm_VA\"" << ","    
                /***************************/
                << "\"sumaRe_VA\"" << "," << "\"sumaIm_VA\"" << ","
                << "\"suma_weRe_VA\"" << "," << "\"suma_weIm_VA\"" << ","    
                /***************************/
                << "\"sumb1Re_VA\"" << "," << "\"sumb1Im_VA\"" << ","
                << "\"sumb1_weRe_VA\"" << "," << "\"sumb1_weIm_VA\"" << ","    
                /***************************/
                << "\"sumb2Re_VA\"" << "," << "\"sumb2Im_VA\"" << ","
                << "\"sumb2_weRe_VA\"" << "," << "\"sumb2_weIm_VA\"" << ","    
                /***************************/
                << "\"sumb3Re_VA\"" << "," << "\"sumb3Im_VA\"" << ","
                << "\"sumb3_weRe_VA\"" << "," << "\"sumb3_weIm_VA\"" << ","    
                /***************************/
                << "\"sumb4Re_VA\"" << "," << "\"sumb4Im_VA\"" << ","
                << "\"sumb4_weRe_VA\"" << "," << "\"sumb4_weIm_VA\"" << ","    
                /***************************/
                << "\"sumb5Re_VA\"" << "," << "\"sumb5Im_VA\"" << ","
                << "\"sumb5_weRe_VA\"" << "," << "\"sumb5_weIm_VA\""
                /***************************/    
                << std::endl;
    
  /* phaseVolume-stream open */
  config.streampc.open(name_files[1], std::ios::out);

  config.streampc << "\"t\"" << "," << "\"T_000Q\"" << ","
                  << "\"Vratio_p1_acc\"" << ","
                  << "\"Vratio_p2_acc\"" << ","
                  << "\"Vratio_p3_acc\"" << ","
                  << "\"Vratio_p4_acc\"" << ","
                  << "\"Vratio_p5_acc\"" << ","
                  << "\"Vratio_p6_acc\"" << ","
                  << "\"Vratio_p7_acc\"" << ","
                  << "\"Vratio_p8_acc\"" << ","
                  << "\"Vratio_p9_acc\"" << ","
                  << "\"Vratio_ps_acc\"" << ","
                  << "\"V_p1_acc\"" << ","
                  << "\"V_p2_acc\"" << ","
                  << "\"V_p3_acc\"" << ","
                  << "\"V_p4_acc\"" << ","
                  << "\"V_p5_acc\"" << ","
                  << "\"V_p6_acc\"" << ","
                  << "\"V_p7_acc\"" << ","
                  << "\"V_p8_acc\"" << ","
                  << "\"V_p9_acc\"" << ","
                  << "\"V_ps_acc\""
		  << std::endl;    
  
    
} // allocate() function ends here

void glsol::fstreams_close() {
  config.stream.close();
  config.streampc.close();
  system("rm -rf pio/*.root pio_Current/*.root 2>/dev/null");
}
  

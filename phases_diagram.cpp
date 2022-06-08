/*
 *
 *  This implement offers datas of the gap parameters, the equlibrium bulk free energies,
 *  
 *  and *zero* magnetic field equlibrium phase diagram of bulk 3He superfluids.
 *
 *  shell script sfdata_plot.sh drives the binary runing and renders the data visualization,
 * 
 *  sfdata.py in shell script calls the matplotlib.
 *
 *
 *  The binary *.app will generate five *csv data files for pressure 0 - 34 bar, T 0 - Tc(p),
 *  the number of points of pressure and temperature are determined by 
 *             
 *           unsigned int num 
 *  
 *  of constructor SFdata( const real_t &, const real_t &, unsigned num ), i.e., line 52.
 *  The dedualt value is set as 400, then one can get 201 pressure points and 401 tempresure points.
 *  
 *  details of sturcture of *csv files will be explained in sfdata.py. 
 *
 *
 *  author: Quang. Zhang (timohyva@github)
 *
 */


#include <iostream>
#include <fstream>
#include <string>
#include <vector>

//#include <cstddef>
//#include <cmath>

#include "matep.hpp"

using real_t = float;
using counter = unsigned int;


/***************************************************************************/
/********                   physics constants                     **********/
/***************************************************************************/
const real_t kb = 1.380649*(std::pow(10.f, -23));


/***************************************************************************/
/********                   SFdata  class                         **********/
/***************************************************************************/

class SFdata {
public:
  SFdata() = default;                                           // default constructor
  SFdata(const real_t &, const real_t &, unsigned num = 400);   // constructor with pressure region,
                                                                // usually the region is [0., 34.] bar. 

  void run();                                                   // driver

private:
  Matep mp;                                             // Matep object, mp

  std::vector<real_t> pressure_list;                    // pressure_list
  std::vector<real_t> unit_list;                        // list of [0. 1.0] for Temperature
  
  std::string name_A;                                   // string for data files
  std::string name_B;
  std::string name_pd;
  
  std::ofstream of_A;                                   // ofstream of_A
  std::ofstream of_B;                                   // ofstream of_B
  std::ofstream of_pd;                                  // ofstream of_pd


  void gaps();                                          // gaps calulator 
  void equlibrium_freeEnergies();                       // equlibrium free energies calculator
  void pd_data();                                       // phase diagram data generator
  
};

/***************************************************************************/
/********                  SFdata constructor                     **********/
/***************************************************************************/  

SFdata::SFdata(const real_t &p_start, const real_t &p_end, unsigned num)
  :pressure_list(num - 199u, 0.)
  ,unit_list(num + 1u, 0.){
  const real_t step_p = (p_end - p_start)/(num - 200u);
  const real_t step_u = 1./num;

  // preseure_list 
  for (auto it_p = pressure_list.begin(); it_p != pressure_list.end(); ++it_p){
    if (it_p == pressure_list.begin()){
      // std::cout << " p = " << *it_p << std::endl;
      continue;
    }else{
      *it_p = *(it_p - 1u) + step_p;
      // std::cout << " p = " << *it_p << std::endl;
    }
  }

  // Tempereture unit_list
  for (auto it_u = unit_list.begin(); it_u != unit_list.end(); ++it_u){
    if (it_u == unit_list.begin()){
      //std::cout << " u_element = " << *it_u << std::endl;
	continue;
      }else{
	*it_u = *(it_u - 1u) + step_u;
      //std::cout << " u_element " << *it_u << std::endl;
      }
  }
}
  

/***************************************************************************/
/********                     gaps function                       **********/
/***************************************************************************/  

void
SFdata::gaps(){

  // ******************   gap_A  ****************** //
  
  name_A = "pd/of_gaps_A.csv";
  of_A.open(name_A);

  for (auto it_p = pressure_list.begin(); it_p != pressure_list.end(); ++it_p){
    of_A << *it_p << ',' << mp.Tcp_mK(*it_p) << ',';
    for (auto it_u =unit_list.begin(); it_u != unit_list.end(); ++it_u){
      auto T = (*it_u) * mp.Tcp_mK(*it_p);
      if (it_u != ((unit_list.end()) - 1u)){
	of_A << mp.gap_A_td(*it_p, T) << ',';
      }else{
	of_A << mp.gap_A_td(*it_p, T);
      }
    }

    of_A << std::endl;
  }

  of_A.close();

  // ******************   gap_B  ****************** //
  
  name_B = "pd/of_gaps_B.csv";
  of_B.open(name_B);

  for (auto it_p = pressure_list.begin(); it_p != pressure_list.end(); ++it_p){
    of_B << *it_p << ',' << mp.Tcp_mK(*it_p) << ',';
    for (auto it_u =unit_list.begin(); it_u != unit_list.end(); ++it_u){
      auto T = (*it_u) * mp.Tcp_mK(*it_p);
      if (it_u != ((unit_list.end()) - 1u)){
	of_B << mp.gap_B_td(*it_p, T) << ',';
      }else{
	of_B << mp.gap_A_td(*it_p, T);
      }
    }

    of_B << std::endl;
  }

  of_B.close();
  
}

/***************************************************************************/
/********                free energies function                   **********/
/***************************************************************************/

void
SFdata::equlibrium_freeEnergies(){

  // *****************  f_A  ******************** //
  
  name_A = "pd/of_f_A.csv";
  of_A.open(name_A);
  
  for (auto it_p = pressure_list.begin(); it_p != pressure_list.end(); ++it_p){
    of_A << *it_p << ',' << mp.Tcp_mK(*it_p) << ',';
    for (auto it_u =unit_list.begin(); it_u != unit_list.end(); ++it_u){
      auto T = (*it_u) * mp.Tcp_mK(*it_p);
      if (it_u != ((unit_list.end()) - 1u)){
	of_A << mp.f_A_td(*it_p, T) << ',';
      }else{
	of_A << mp.f_A_td(*it_p, T);
      }
    }

    of_A << std::endl;
  }

  of_A.close();

  // *****************  f_B  ******************** //
  
  name_B = "pd/of_f_B.csv";
  of_B.open(name_B);
  
  for (auto it_p = pressure_list.begin(); it_p != pressure_list.end(); ++it_p){
    of_B << *it_p << ',' << mp.Tcp_mK(*it_p) << ',';
    for (auto it_u =unit_list.begin(); it_u != unit_list.end(); ++it_u){
      auto T = (*it_u) * mp.Tcp_mK(*it_p);
      if (it_u != ((unit_list.end()) - 1u)){
	of_B << mp.f_B_td(*it_p, T) << ',';
      }else{
	of_B << mp.f_B_td(*it_p, T);
      }
    }

    of_B << std::endl;
  }

  of_B.close();
}

/***************************************************************************/
/********                     pd_data function                     *********/
/***************************************************************************/

void
SFdata::pd_data(){
  name_pd = "pd/phase_diagram_data.csv";
  of_pd.open(name_pd);

  for (auto it_p = pressure_list.begin(); it_p != pressure_list.end(); ++it_p){
    if (*it_p >= 21.22){
     of_pd << *it_p << ','
           << mp.tAB_RWS(*it_p)*mp.Tcp_mK(*it_p) << ','
	   << mp.Tcp_mK(*it_p);
    }else
      of_pd << *it_p << ','
           << 0.0f << ','
	   << mp.Tcp_mK(*it_p);
    
    of_pd << std::endl;
  }

  of_pd.close();
}


/***************************************************************************/
/********                     driver function                     **********/
/***************************************************************************/  

void
SFdata::run(){
  gaps();
  equlibrium_freeEnergies();
  pd_data();
}


/***************************************************************************/
/********                           main                          **********/
/***************************************************************************/  


int main(){

  const real_t p_0 = 0.f, p_34 = 34.f;    // bar
  
  SFdata sf_data(p_0, p_34);

  sf_data.run();

  return 0;

}  

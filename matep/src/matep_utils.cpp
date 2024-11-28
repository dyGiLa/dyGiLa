/*
 * This is the *.cpp file of the Strong Coupling Correction Object (SCCO)
 * or Material Parameters Object (Matep).
 * 
 * Member functions are declared at Matep.hpp, and then defined at here.
 *
 * author: Quang. Zhang (timohyva@github)
 *
 */ 


#include <iostream>
#include <cstddef>
#include <cmath>
#include <vector>

//#include "matep_namespace_utils.hpp"
#include "matep.hpp"

namespace matep {


//**********************************************************************
//***              public member: Levi-Civita symbol                 ***
//**********************************************************************

real_t
Matep::epsilon(int al, int be, int ga)
{
  if (
      (al == 0 && be == 1 && ga == 2)   // 123
      ||(al == 1 && be == 2 && ga == 0) // 231
      ||(al == 2 && be == 0 && ga == 1) // 312
     )
    {return 1.0;}
  else if
     (
      (al == 2 && be == 1 && ga == 0)   // 321
      ||(al == 0 && be == 2 && ga == 1) // 132
      ||(al == 1 && be == 0 && ga == 2) // 213
     )
    {return -1.0;}
  else if ((al == be) || (al == ga) || (be == ga))
    {return 0.0;}
  else
    {return 0.0;}
}  
  

//**********************************************************************
//***       private method :  linear intepolation function           ***
//**********************************************************************

real_t
Matep::lininterp(const real_t *cX_arr, real_t p){
  float pk, pk1, fp;
  size_t k, k1;

  if ((p >= 0.0) && (p < 2.0)) { pk = 0.0f; k = 0; }

  if ((p >= 2.0) && (p < 4.0)) { pk = 2.0f; k = 1; }
 
  if ((p >= 4.0) && (p < 6.0)) { pk = 4.0f; k = 2; }

  if ((p >= 6.0) && (p < 8.0)) { pk = 6.0f; k = 3; }

  if ((p >= 8.0) && (p < 10.0)) { pk = 8.0f; k = 4; }

  if ((p >= 10.0) && (p < 12.0)) { pk = 10.0f; k = 5; }

  if ((p >= 12.0) && (p < 14.0)) { pk = 12.0f; k = 6; }

  if ((p >= 14.0) && (p < 16.0)) { pk = 14.0f; k = 7; }

  if ((p >= 16.0) && (p < 18.0)) { pk = 16.0f; k = 8; }

  if ((p >= 18.0) && (p < 20.0)) { pk = 18.0f; k = 9; }

  if ((p >= 20.0) && (p < 22.0)) { pk = 20.0f; k = 10; }

  if ((p >= 22.0) && (p < 24.0)) { pk = 22.0f; k = 11; }

  if ((p >= 24.0) && (p < 26.0)) { pk = 24.0f; k = 12; }

  if ((p >= 26.0) && (p < 28.0)) { pk = 26.0f; k = 13; }

  if ((p >= 28.0) && (p < 30.0)) { pk = 28.0f; k = 14; }

  if ((p >= 30.0) && (p < 32.0)) { pk = 30.0f; k = 15; }

  if ((p >= 32.0) && (p < 34.0)) { pk = 32.0f; k = 16; }

  fp = ((cX_arr[k+1]-cX_arr[k])/2.0)*(p-pk)+cX_arr[k];
  return fp; 
}

} // namespace matep block ends here  


//////////////////////
//R
//1D Field under Polar co

#ifndef GEFICA_R_H
#define GEFICA_R_H

#include "X.h"

namespace GEFICA { 
   class R;
}

class GEFICA::R : public X 
{
   public:
      R(int nx=101): X(nx){};
      ClassDef(R,1);

   protected:

      virtual void SOR2(int idx,bool elec); // 2nd-order Runge-Kutta Successive Over-Relaxation
      virtual void SOR4(int idx); // 4th-order Runge-Kutta Successive Over-Relaxation
};

#endif

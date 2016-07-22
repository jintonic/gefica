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
      virtual void RK2(int idx); // 2nd-order Runge-Kutta Successive Over-Relaxation
      virtual void RK4(int idx); // 4th-order Runge-Kutta Successive Over-Relaxation
};

#endif

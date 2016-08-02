#ifndef GEFICA_Rho_H
#define GEFICA_Rho_H

#include "X.h"

namespace GEFICA { 
   class Rho;
}

class GEFICA::Rho : public X 
{
   public:
      Rho(int nx=101): X(nx){};
      ClassDef(Rho,1);

   protected:

      virtual void RK2(int idx,bool elec); // 2nd-order Runge-Kutta Successive Over-Relaxation
      virtual void RK4(int idx); // 4th-order Runge-Kutta Successive Over-Relaxation
};

#endif

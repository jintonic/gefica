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
      virtual double GetPotential(double rho){return GetData(rho,1);};
      virtual double GetE1(double rho){return GetData(rho,2);};
      virtual double GetImpurity(double rho){return GetData(rho,0);};
   protected:

      virtual void SOR2(int idx,bool elec); // 2nd-order Runge-Kutta Successive Over-Relaxation
      virtual void SOR4(int idx); // 4th-order Runge-Kutta Successive Over-Relaxation
};

#endif

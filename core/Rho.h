#ifndef GeFiCa_Rho_H
#define GeFiCa_Rho_H

#include "X.h"

namespace GeFiCa { 
   class Rho;
}

class GeFiCa::Rho : public X 
{
   public:
      /**
	  * Rho is a constructor, if given a number, no input is needed
	  */
      Rho(int nx=101): X(nx){};
	  /**
	  * This defines the class R for the cint dictionary.
	  */
      ClassDef(Rho,1);
      virtual double GetPotential(double rho){return GetData(rho,1);};
      virtual double GetE1(double rho){return GetData(rho,2);};
      virtual double GetImpurity(double rho){return GetData(rho,0);};
   protected:

      virtual void SOR2(int idx,bool elec); // 2nd-order Runge-Kutta Successive Over-Relaxation
      virtual void SOR4(int idx); ///< 4th-order Runge-Kutta Successive Over-Relaxation
};
#endif


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
      Rho(int nx=101): X(nx){t=2,d=1;};
	  /**
	  * This defines the class R for the cint dictionary.
	  */
      ClassDef(Rho,1);
      virtual double GetPotential(double rho){return GetData(rho,kPotential);};
      virtual double GetE1(double rho,double phi,double z){return GetData(rho,kE1);};
      virtual double GetImpurity(double rho){return GetData(rho,kImpurity);};
   protected:

      virtual void SOR2(int idx,bool elec); // 2nd-order Runge-Kutta Successive Over-Relaxation
      virtual void SOR4(int idx); ///< 4th-order Runge-Kutta Successive Over-Relaxation
};
#endif


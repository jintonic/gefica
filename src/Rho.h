#ifndef GeFiCa_Rho_H
#define GeFiCa_Rho_H

#include "X.h"

namespace GeFiCa { class Rho; }

/**
 * 1D cylindrical coordinate.
 */
class GeFiCa::Rho : public X 
{
   public:
      /**
       * Default constructor.
       */
      Rho(int n=101, const char *name="rho",
            const char *title="1D cylindrical coordinates")
         : X(n, name, title) {};

      virtual double GetPotential(double rho)
      {return GetData(rho,kPotential);};
      virtual double GetE1(double rho,double phi,double z)
      {return GetData(rho,kE1);};
      virtual double GetImpurity(double rho)
      {return GetData(rho,kImpurity);};

      ClassDef(Rho,1);

   protected:
      virtual void DoSOR2(int idx);
      virtual void DoSOR4(int idx);
};
#endif


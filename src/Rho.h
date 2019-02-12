#ifndef GeFiCa_Rho_H
#define GeFiCa_Rho_H

#include "X.h"

namespace GeFiCa { 
   class Rho;
}

class GeFiCa::Rho : public X 
{
   public:
      Rho(int nx=101): X(nx) {SetName("Rho"); SetTitle("Rho"); }

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


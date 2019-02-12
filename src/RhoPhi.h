#ifndef GeFiCa_RHOPHI_H
#define GeFiCa_RHOPHI_H

#include "XY.h"

namespace GeFiCa { class RhoPhi; }

class GeFiCa::RhoPhi : public GeFiCa::XY
{
   public:
      RhoPhi(int nRho=101, int nPhi=101): XY(nRho,nPhi) {}; 

      virtual double GetPotential(double rho,double phi)
      {return GetData(rho,phi,kPotential);};
      virtual double GetE1(double rho,double phi,double z)
      {return GetData(rho,phi,kE1);};
      virtual double GetE2(double rho,double phi,double z)
      {return GetData(rho,phi,kE2);};
      virtual double GetImpurity(double rho,double phi)
      {return GetData(rho,phi,kImpurity);};

      ClassDef(RhoPhi,1);
   protected:
      void DoSOR2(int idx);
      virtual double GetData(double tarx,double tary,EOutput output);
};
#endif


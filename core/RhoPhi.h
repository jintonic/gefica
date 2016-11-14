#ifndef GeFiCa_RHOPHI_H
#define GeFiCa_RHOPHI_H

#include "XY.h"

namespace GeFiCa { class RhoPhi; }

class GeFiCa::RhoPhi : public GeFiCa::XY
{
   public:
      RhoPhi(unsigned short r=101, unsigned short phi=101): XY(r,phi) {} 
      virtual ~RhoPhi() {};

      virtual double GetPotential(double rho,double phi){return GetData(rho,phi,1);};
      virtual double GetE1(double rho,double phi){return GetData(rho,phi,2);};
      virtual double GetE2(double rho,double phi){return GetData(rho,phi,3);};
      virtual double GetImpurity(double rho,double phi){return GetData(rho,phi,0);};
   protected:
      virtual double GetData(double tarx,double tary,int thing);
      void SOR2(int idx,bool elec);

      ClassDef(RhoPhi,1);
};
#endif


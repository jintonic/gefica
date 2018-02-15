#ifndef GeFiCa_RHOPHI_H
#define GeFiCa_RHOPHI_H

#include "XY.h"

namespace GeFiCa { class RhoPhi; }

class GeFiCa::RhoPhi : public GeFiCa::XY
{
   public:
      RhoPhi(unsigned short r=101, unsigned short phi=101): XY(r,phi) {t=2,d=2;}; 
      virtual ~RhoPhi() {};

      virtual double GetPotential(double rho,double phi){return GetData(rho,phi,kPotential);};
      virtual double GetE1(double rho,double phi,double z){return GetData(rho,phi,kE1);};
      virtual double GetE2(double rho,double phi,double z){return GetData(rho,phi,kE2);};
      virtual double GetImpurity(double rho,double phi){return GetData(rho,phi,kImpurity);};
   protected:
      virtual double GetData(double tarx,double tary,EOutput output);
      void SOR2(int idx,bool elec);

      ClassDef(RhoPhi,1);
};
#endif


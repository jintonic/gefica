#ifndef GeFiCa_RHOPHIZ_H
#define GeFiCa_RHOPHIZ_H

#include "XYZ.h"
class TF3;

namespace GeFiCa { class RhoPhiZ; }

class GeFiCa::RhoPhiZ : public GeFiCa::XYZ
{
   public:
      RhoPhiZ(unsigned short n1, unsigned short n2,unsigned short n3): XYZ(n1,n2,n3) {};
      virtual ~RhoPhiZ(){};
      virtual double GetPotential(double rho,double phi,double z){return GetData(rho,phi,z,1);};
      virtual double GetE1(double rho,double phi,double z){return GetData(rho,phi,z,2);};
      virtual double GetE2(double rho,double phi,double z){return GetData(rho,phi,z,3);};
      virtual double GetE3(double rho,double phi,double z){return GetData(rho,phi,z,4);};
      virtual double GetImpurity(double rho,double phi,double z){return GetData(rho,phi,z,0);};
   protected:
      virtual void SOR2(int idx,bool elec); 

      virtual double GetData(double tarx,double tary,double tarz,int thing);

      ClassDef(RhoPhiZ,1);

};

#endif

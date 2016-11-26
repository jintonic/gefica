#ifndef GeFiCa_RTHETAPHI_H
#define GeFiCa_RTHETAPHI_H

#include "XYZ.h"

namespace GeFiCa { class RThetaPhi; }

class GeFiCa::RThetaPhi : public GeFiCa::XYZ
{

   public:
      RThetaPhi(unsigned short n_r=0, unsigned short n_theta=0,unsigned short n_phi=0): 
         XYZ(n_r, n_theta, n_phi*2) {t=3;d=3;}
      virtual ~RThetaPhi() {};
      virtual double GetPotential(double r,double theta,double phi){return GetData(r,theta,phi,1);};
      virtual double GetE1(double r,double theta,double phi){return GetData(r,theta,phi,2);};
      virtual double GetE2(double r,double theta,double phi){return GetData(r,theta,phi,3);};
      virtual double GetE3(double r,double theta,double phi){return GetData(r,theta,phi,4);};
      virtual double GetImpurity(double r,double theta,double phi){return GetData(r,theta,phi,0);};
   protected:
      virtual void SOR2(int idx,bool elec); 

      virtual double GetData(double tarx,double tary,double tarz,int thing);
      ClassDef(RThetaPhi,1);

};
#endif


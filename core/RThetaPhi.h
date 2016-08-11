#ifndef GEFICA_RTHETAPHI_H
#define GEFICA_RTHETAPHI_H

#include "XYZ.h"

namespace GEFICA { class RThetaPhi; }

class GEFICA::RThetaPhi : public GEFICA::XYZ
{

   public:
      RThetaPhi(unsigned short n_r=0, unsigned short n_theta=0,unsigned short n_phi=0): 
         XYZ(n_r, n_theta, n_phi*2) {}
      virtual ~RThetaPhi() {};

      virtual void SOR2(int idx,bool elec); 

      virtual double GetData(double tarx,double tary,double tarz,int thing);
      ClassDef(RThetaPhi,1);

};

#endif

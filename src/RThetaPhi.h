#ifndef GeFiCa_RTHETAPHI_H
#define GeFiCa_RTHETAPHI_H

#include "XYZ.h"

namespace GeFiCa { class RThetaPhi; }

class GeFiCa::RThetaPhi : public GeFiCa::XYZ
{
   public:
      RThetaPhi(int n_r=0, int n_theta=0, int n_phi=0): 
         XYZ(n_r, n_theta, n_phi*2) {};

      ClassDef(RThetaPhi,1);

   protected:
      virtual void DoSOR2(int idx);

      virtual double GetData(double tarx,double tary,double tarz, EOutput output);
};
#endif


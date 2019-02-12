#ifndef GeFiCa_RHOPHIZ_H
#define GeFiCa_RHOPHIZ_H

#include "XYZ.h"
class TF3;

namespace GeFiCa { class RhoPhiZ; }

class GeFiCa::RhoPhiZ : public GeFiCa::XYZ
{
   public:
      RhoPhiZ(int n1, int n2, int n3): XYZ(n1,n2,n3) {};

      ClassDef(RhoPhiZ,1);

   protected:
      virtual void DoSOR2(int idx); 

      virtual double GetData(double tarx,double tary,double tarz, EOutput output);
};

#endif

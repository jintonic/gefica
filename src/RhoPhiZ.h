#ifndef GeFiCa_RHOPHIZ_H
#define GeFiCa_RHOPHIZ_H

#include "XYZ.h"

namespace GeFiCa { class RhoPhiZ; }

/**
 * 3D cylindrical coordinates.
 */
class GeFiCa::RhoPhiZ : public GeFiCa::XYZ
{
   public:
      RhoPhiZ(int fN1, int fN2, int fN3): XYZ(fN1,fN2,fN3) {};

      ClassDef(RhoPhiZ,1);

   protected:
      virtual void DoSOR2(int idx); 

      virtual double GetData(double tarx,double tary,double tarz, EOutput output);
};

#endif

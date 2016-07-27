#ifndef GEFICA_TRUEC_H
#define GEFICA_TRUEC_H

#include "Cylindrical.h"

namespace GEFICA { class TrueC; }

class GEFICA::TrueC : public GEFICA::Cylindrical
{
   public :
      double r1,r0;

      TrueC(int r,int O,int z) : Cylindrical(r, O, z) {};
      void SetVoltage(double v0,double v1); 
      void Create(double r0,double r1);

      ClassDef(TrueC,1);
};

#endif

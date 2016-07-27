#ifndef GEFICA_TRUEC2D_H
#define GEFICA_TRUEC2D_H

#include "Polar.h"

namespace GEFICA { class TrueC2d; }

class GEFICA::TrueC2d : public GEFICA::Polar
{
   public :
      double r1,r0;

      TrueC2d(int r,int O) : Polar(r, O) {};
      void SetVoltage(double v0,double v1); 
      void Create(double r0,double r1);

      ClassDef(TrueC2d,1);
};

#endif

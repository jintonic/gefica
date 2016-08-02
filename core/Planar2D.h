#ifndef GEFICA_PLANAR2D_H
#define GEFICA_PLANAR2D_H

#include "XY.h"

namespace GEFICA { class Planar2D; }

class GEFICA::Planar2D : public GEFICA::XY
{
   public :
      double Thickness, Width;

   public :
      Planar2D(int ix,int iy) : XY(ix,iy) {};
      void SetVoltage(double voltage); 

      ClassDef(Planar2D, 1);
};

#endif

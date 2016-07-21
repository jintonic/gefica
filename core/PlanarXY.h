#ifndef GEFICA_PLANARXY_H
#define GEFICA_PLANARXY_H

#include "XY.h"

namespace GEFICA { class PlanarXY; }

class GEFICA::PlanarXY : public GEFICA::XY
{
   public :
      double Thickness, Width;

   public :
      PlanarXY(int ix,int iy) : XY(ix,iy) {};
      void SetVoltage(double voltage); 

      ClassDef(PlanarXY, 1);
};

#endif

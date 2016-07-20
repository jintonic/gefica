#ifndef GEFICA_PLANARXY_H
#define GEFICA_PLANARXY_H

#include "XY.h"

namespace GEFICA { class PlanarXY; }

class GEFICA::PlanarXY : public GEFICA::XY
{
   public :
      bool isx;

   public :
      PlanarXY(int ix,int iy) : XY(ix,iy) {};
      void SetVoltage(double voltage); 

      double GetWidth() { return fC2[n]; }
};

#endif

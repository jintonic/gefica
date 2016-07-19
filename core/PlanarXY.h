#ifndef PLANARXY_H
#define PLANARXY_H

#include "XY.h"

namespace GEFICA { class PlanarXY; }

class GEFICA::PlanarXY : public GEFICA::XY
{
  public :
    bool isx;
    PlanarXY(int ix,int iy) : XY(ix,iy) {};
    void SetVoltage(double voltage); 
    ClassDef(PlanarXY,1);
};

#endif

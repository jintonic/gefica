#ifndef PLANAR2D_H
#define PLANAR2D_H

#include "Field2D.h"

namespace GEFICA { class Planar2D; }

class GEFICA::Planar2D : public GEFICA::Field2D
{
  public :
    bool isx;
    Planar2D(int ix,int iy) : Field2D(ix,iy) {};
    void SetVoltage(double voltage); 
    ClassDef(Planar2D,1);
};

#endif

#ifndef NSHAPE_H
#define NSHAPE_H

#include "Field2D.h"

namespace GEFICA { class nshape; }

class GEFICA::nshape : public GEFICA::Field2D
{
  public :
    nshape(int ix,int iy) : Field2D(ix,iy) {};
    void SetVoltage(double voltage,double topbegin,double topend); 
    ClassDef(nshape,1);
};

#endif

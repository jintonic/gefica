#ifndef PLANAR1D_H
#define PLANAR1D_H

#include "Field.h"

namespace GEFICA { class Planar1D; }

class GEFICA::Planar1D : public GEFICA::Field
{
  public :
    Planar1D(int ix) : Field(ix) {};
    void SetVoltage(double voltage); 

    ClassDef(Planar1D,1);
};

#endif

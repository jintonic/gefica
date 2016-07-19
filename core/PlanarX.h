#ifndef PLANARX_H
#define PLANARX_H

#include "X.h"

namespace GEFICA { class PlanarX; }

class GEFICA::PlanarX : public GEFICA::X
{
  public :
    PlanarX(int ix) : X(ix) {};
    void SetVoltage(double voltage); 

    ClassDef(PlanarX,1);
};

#endif

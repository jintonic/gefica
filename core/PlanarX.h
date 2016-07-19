#ifndef GEFICA_PLANARX_H
#define GEFICA_PLANARX_H

#include "X.h"

namespace GEFICA { class PlanarX; }

class GEFICA::PlanarX : public GEFICA::X
{
   public :
      PlanarX(int nx=11) : X(nx) {};
      void SetVoltage(double voltage); 
};

#endif

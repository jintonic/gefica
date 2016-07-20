#ifndef GEFICA_PLANARX_H
#define GEFICA_PLANARX_H

#include "X.h"

namespace GEFICA { class PlanarX; }

class GEFICA::PlanarX : public GEFICA::X
{
   public :
      PlanarX(int nx=11) : X(nx) {};
      double GetThickness() { return fC1[n]; }
      void SetVoltage(double anode_voltage, double cathode_voltage); 
};

#endif

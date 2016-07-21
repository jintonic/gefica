#ifndef GEFICA_PLANARX_H
#define GEFICA_PLANARX_H

#include "X.h"

namespace GEFICA { class PlanarX; }

class GEFICA::PlanarX : public GEFICA::X
{
   public :
      double Thickness; // Thickness of the planar detector

   public :
      PlanarX(int nx=101) : X(nx) {};
      void SetVoltage(double anode_voltage, double cathode_voltage); 
      virtual bool Analyic();
      
      ClassDef(PlanarX, 1);
};

#endif

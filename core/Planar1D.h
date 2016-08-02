#ifndef GEFICA_PLANAR1D_H
#define GEFICA_PLANAR1D_H

#include "X.h"

namespace GEFICA { class Planar1D; }

class GEFICA::Planar1D : public GEFICA::X
{
   public :
      double Thickness; // Thickness of the planar detector

   public :
      Planar1D(int nx=101) : X(nx) {};
      void SetVoltage(double anode_voltage, double cathode_voltage); 
      bool Analyic();
      
      ClassDef(Planar1D, 1);
};

#endif

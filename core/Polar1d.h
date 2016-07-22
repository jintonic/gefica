#ifndef GEFICA_POLAR1D_H
#define GEFICA_POLAR1D_H

#include "R.h"

namespace GEFICA { class Polar1d; }

class GEFICA::Polar1d : public GEFICA::R
{
   public :
      double Thickness; // Thickness of the planar detector

   public :
      Polar1d(int nx=101) : R(nx) {};
      void SetVoltage(double anode_voltage, double cathode_voltage); 
//      bool Analyic();
      
      ClassDef(Polar1d, 1);
};

#endif

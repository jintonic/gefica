#ifndef GEFICA_TRUECOAXIAL1D_H
#define GEFICA_TRUECOAXIAL1D_H

#include "R.h"

namespace GEFICA { class TrueCoaxial1D; }

class GEFICA::TrueCoaxial1D : public GEFICA::R
{
   public :
      double Thickness; // Thickness of the planar detector

   public :
      TrueCoaxial1D(int nx=101) : R(nx) {};
      void SetVoltage(double anode_voltage, double cathode_voltage); 
      bool Analyic();
      void Create(double r0,double r1);
      ClassDef(TrueCoaxial1D, 1);
};

#endif

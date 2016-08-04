#ifndef GEFICA_SPHERE1D_H
#define GEFICA_SPHERE1D_H

#include "R.h"

namespace GEFICA { class Sphere1D; }

class GEFICA::Sphere1D : public GEFICA::R
{
   public :
      double InnerRadius; // Inner radius of the detector
      double OuterRadius; // Outer radius of the detector

   public :
      Sphere1D(int nx=101) : R(nx) {};
      void SetVoltage(double anode_voltage, double cathode_voltage); 
      void CreateGridWithFixedStepLength();

      ClassDef(Sphere1D, 1);

   protected:
      bool Analyic();
};

#endif

#ifndef GEFICA_HEMISPHERIC_H
#define GEFICA_HEMISPHERIC_H

#include "Rho.h"

namespace GEFICA { class HemiSpheric; }

class GEFICA::HemiSpheric : public GEFICA::Rho
{
   public :
      double Thickness; // Thickness of the planar detector

   public :
      HemiSpheric(int nx=101) : Rho(nx) {};
      void SetVoltage(double anode_voltage, double cathode_voltage); 
      bool Analyic();
      void Create(double r0,double r1);
      ClassDef(HemiSpheric, 1);
};

#endif

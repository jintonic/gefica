#ifndef GEFICA_SPH_H
#define GEFICA_SPH_H

#include "Rho.h"

namespace GEFICA { class Sph; }

class GEFICA::Sph : public GEFICA::Rho
{
   public :
      double Thickness; // Thickness of the planar detector

   public :
      Sph(int nx=101) : Rho(nx) {};
      void SetVoltage(double anode_voltage, double cathode_voltage); 
      bool Analyic();
      void Create(double r0,double r1);
      ClassDef(Sph, 1);
};

#endif

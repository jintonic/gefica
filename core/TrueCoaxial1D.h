#ifndef GEFICA_TRUECOAXIAL1D_H
#define GEFICA_TRUECOAXIAL1D_H

#include "Rho.h"

namespace GEFICA { class TrueCoaxial1D; }

class GEFICA::TrueCoaxial1D : public GEFICA::Rho
{
   public :
      double OuterRadius; // Outer radius of the detector
      double InnerRadius; // Inner radius of the detector
      double cathode_voltage,annode_voltage;

   public :
      TrueCoaxial1D(int nx=101) : Rho(nx),OuterRadius(10),InnerRadius(0), cathode_voltage(2000),annode_voltage(0) {};
      void initialize();
      bool CalculateField(EMethod method=kSOR2);
      ClassDef(TrueCoaxial1D, 1);

   protected:
      bool Analytic();
};

#endif

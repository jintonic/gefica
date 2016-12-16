#ifndef GeFiCa_TRUECOAXIAL1D_H
#define GeFiCa_TRUECOAXIAL1D_H

#include "Rho.h"

namespace GeFiCa { class TrueCoaxial1D; }

class GeFiCa::TrueCoaxial1D : public GeFiCa::Rho
{
   public :
      double OuterRadius; // Outer radius of the detector
      double InnerRadius; // Inner radius of the detector

   public :
      TrueCoaxial1D(int nx=101) : Rho(nx),OuterRadius(10),InnerRadius(0) {};
      void Initialize();
      bool CalculateField(EMethod method=kSOR2);
      ClassDef(TrueCoaxial1D, 1);

   protected:
      bool Analytic();
};
#endif


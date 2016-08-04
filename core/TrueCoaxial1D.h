#ifndef GEFICA_TRUECOAXIAL1D_H
#define GEFICA_TRUECOAXIAL1D_H

#include "Rho.h"

namespace GEFICA { class TrueCoaxial1D; }

class GEFICA::TrueCoaxial1D : public GEFICA::Rho
{
   public :
      double InnerRadius; // Inner radius of the detector
      double OuterRadius; // Outer radius of the detector

   public :
      TrueCoaxial1D(int nx=101) : Rho(nx) {};

      ClassDef(TrueCoaxial1D, 1);

   protected:
      bool Analyic();
};

#endif

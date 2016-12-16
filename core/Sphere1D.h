#ifndef GeFiCa_SPHERE1D_H
#define GeFiCa_SPHERE1D_H

#include "R.h"

namespace GeFiCa { class Sphere1D; }

class GeFiCa::Sphere1D : public GeFiCa::R
{
  public:
      double innerR,outterR; // boundary of the planar detector

   public :
      Sphere1D(int nx=101) : R(nx),innerR(0.3),outterR(3) {};
      ClassDef(Sphere1D, 1);
      void Initialize();      
      bool CalculateField(EMethod method=kSOR2);
   protected:
      bool Analytic();
};

#endif

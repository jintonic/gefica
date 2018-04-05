#ifndef GeFiCa_SPHERE1D_H
#define GeFiCa_SPHERE1D_H

#include "R.h"

namespace GeFiCa { class Sphere1D; }

class GeFiCa::Sphere1D : public GeFiCa::R
{
  public:
      double InnerRadius,OuterRadius; // boundary of the planar detector

   public :
      Sphere1D(int nx=101) : R(nx),InnerRadius(0.3),OuterRadius(3) {};
      ClassDef(Sphere1D, 1);
      void Initialize();      
      bool CalculatePotential(EMethod method=kSOR2);
   protected:
      /**
       * Analytic calculation of 1D field in spheric coordinates.
       *
       * According to 
       * https://www.wolframalpha.com/input/?i=f%28x%29%27%27+%2B2%2Fx*f%28x%29%27%2Ba%3D0
       * In case of fixed impurity, potential(x) = -rho/6/epsilon*r^2 + c1/r +
       * c2 with boundary conditions:
       *
       * - potential(InnerRadius) = V0,
       * - potential(OuterRadius) = V1,
       *
       * So, 
       * - c1 = 
       * - c2 =
       */
      bool Analytic();
};

#endif

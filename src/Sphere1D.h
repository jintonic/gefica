#ifndef GeFiCa_SPHERE1D_H
#define GeFiCa_SPHERE1D_H

#include "R.h"

namespace GeFiCa { class Sphere1D; }

/**
 * Grid setup for 1D spherical detectors.
 */
class GeFiCa::Sphere1D : public GeFiCa::R
{
  public:
      double InnerRadius; ///< inner radius
      double OuterRadius; ///< outer radius

      /**
       * Default constructor.
       */
      Sphere1D(int n=101, const char *name="s1d",
            const char *title="1D spherical detector");

      virtual void Initialize();      

      ClassDef(Sphere1D, 1);

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

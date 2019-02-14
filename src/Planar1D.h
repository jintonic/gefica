#ifndef GeFiCa_PLANAR1D_H
#define GeFiCa_PLANAR1D_H

#include "X.h"

namespace GeFiCa { class Planar1D; }
/**
 * Configuration for 1D planar detectors.
 */
class GeFiCa::Planar1D : public GeFiCa::X
{
   public :
      double Thickness; ///< thickness of planar detector

      Planar1D(int n=101); ///< default constructor

      void Initialize(); ///< initialize data members

      bool CalculatePotential(EMethod method=kSOR2);

      ClassDef(Planar1D, 1);

   protected:
      /**
       * Analytic calculation of 1D field with fixed impurity concentration.
       *
       * In case of fixed impurity, rho, the solution of Poisson's Equation
       *
       *     d^2 p / dx^2 = -rho/epsilon
       *
       * with boundary conditions:
       *
       * - potential(0) = V0,
       * - potential(d) = V1,
       *
       * where, d = UpperBound - LowerBound, is
       *
       *     potential(x) = a x^2 + b x + c 
       *
       * where, 
       *
       * - a = - rho/2/epsilon
       * - b = (V1-V0 - ad^2)/d
       * - c = V0
       */
      bool Analytic();
};
#endif


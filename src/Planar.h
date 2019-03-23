#ifndef GeFiCa_Planar
#define GeFiCa_Planar

#include <TNamed.h>

#include "X.h"
#include "XY.h"
#include "Detector.h"

namespace GeFiCa { class Planar; }
/**
 * Grid setup for planar detectors.
 */
class GeFiCa::Planar : public GeFiCa::Detector, public TNamed
{
   public :
      double Width; ///< width of planar detector

      Planar() : Detector(), TNamed("planar", "planar detector") {};
      /**
       * Analytic calculation of  field with fixed impurity concentration.
       *
       * In case of fixed impurity, rho, the solution of Poisson's Equation
       *
       *     d^2 p / dx^2 = -rho/epsilon
       *
       * with boundary conditions:
       *
       * - potential(0) = Bias[0],
       * - potential(d) = Bias[1],
       *
       * where, d = UpperBound - LowerBound, is
       *
       *     potential(x) = a x^2 + b x + c 
       *
       * where, 
       *
       * - a = - rho/2/epsilon
       * - b = (Bias[1]-Bias[0] - ad^2)/d
       * - c = Bias[0]
       */
      void FillGridWithAnalyticResult();

      void Configure(Grid& grid)
      { if (N2==0) ConfigureX(X& grid); else ConfigureXY(XY& grid); }

      ClassDef(Planar, 1);

   protected:
      void ConfigureX(X& grid) {};
      void ConfigureXY(XY& grid) {};
};
#endif


#ifndef GeFiCa_Hemispherical
#define GeFiCa_Hemispherical

#include <TNamed.h>

#include "R.h"
#include "RTheta.h"
#include "Detector.h"

namespace GeFiCa { class Hemispherical; }
/**
 * Grid setup for spherical detectors.
 */
class GeFiCa::Hemispherical : public GeFiCa::Detector, public TNamed
{
  public:
      double PointContactR; ///< radius of point contact
      double PointContactH; ///< height of point contact

      Hemispherical() : Detector(), TNamed("hsd", "hemispherical detector") {};
      /**
       * Analytic calculation of 1D field in spheric coordinates.
       * According to 
       * https://www.wolframalpha.com/input/?i=f%28x%29%27%27+%2B2%2Fx*f%28x%29%27%2Ba%3D0
       * In case of fixed impurity,
       *    potential(x) = -rho/6/epsilon*r^2 + c1/r + c2
       * with boundary conditions:
       * - potential(InnerR) = Bias[0],
       * - potential(OuterR) = Bias[1],
       */
      void FillGridWithAnalyticResult();

      void Configure(Grid& grid)
      { if (N2==0) ConfigureR(RX& grid); else Configure2D(RTheta& grid); }

      ClassDef(Hemispherical, 1);

   protected:
      void ConfigureR(RX& grid) {};
      void Configure2D(RTheta& grid) {};
};

#endif

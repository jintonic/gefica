#ifndef GeFiCa_X
#define GeFiCa_X
#include <TNamed.h>
#include "Grid.h"
namespace GeFiCa { class X; }
/**
 * 1D Cartesian coordinate.
 */
class GeFiCa::X : public GeFiCa::Grid, public TNamed
{
   public:
      X(size_t n1=101) : Grid(n1), TNamed("x","1D coordinate") {};
      void SetBoundaryCondition(Detector *detector);
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

      ClassDef(X,1);
   protected:
      void OverRelaxAt(size_t idx); 
};
#endif


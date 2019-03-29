#ifndef GeFiCa_X
#define GeFiCa_X

#include <TNamed.h>
#include "Grid.h"
namespace GeFiCa { class X; }
/**
 * 1D Cartesian coordinate.
 */
class GeFiCa::X : public Grid, public TNamed
{
   public:
      X(size_t n1=101) : Grid(n1), TNamed("x","1D Cartesian coordinate") {};

      void GetBoundaryConditionFrom(Detector &detector);
      void SolveAnalytically();

      ClassDef(X,1);
   protected:
      void OverRelaxAt(size_t idx); 
      void CalculateE();
};
#endif

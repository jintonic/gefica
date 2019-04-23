#ifndef GeFiCa_R
#define GeFiCa_R

#include <TNamed.h>
#include "Grid.h"
namespace GeFiCa { class R; }
/**
 * 1D spherical coordinate.
 */
class GeFiCa::R : public Grid, public TNamed
{
   public:
      R(size_t n1=101) : Grid(n1),
      TNamed("r", "1D spherical coordinate") {};

      void GetBoundaryConditionFrom(Detector &detector);
      void SolveAnalytically();
      double GetC();

      ClassDef(R,1);
   protected:
      virtual void OverRelaxAt(size_t idx);
};
#endif


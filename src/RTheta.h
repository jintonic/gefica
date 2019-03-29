#ifndef GeFiCa_RTheta
#define GeFiCa_RTheta

#include <TNamed.h>
#include "Grid.h"
namespace GeFiCa { class RTheta; }
/**
 * 2D spherical coordinates.
 */
class GeFiCa::RTheta : public Grid, public TNamed
{
   public:
      RTheta(size_t n1=101, size_t n2=181) : Grid(n1, n2),
      TNamed("rt", "2D spherical coordinate") {};

      void GetBoundaryConditionFrom(Detector &detector);

      ClassDef(RTheta,1);
   protected:
      virtual void OverRelaxAt(size_t idx);
      void CalculateE();
};
#endif


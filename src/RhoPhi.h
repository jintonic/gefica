#ifndef GeFiCa_RhoPhi
#define GeFiCa_RhoPhi

#include <TNamed.h>
#include "Grid.h"
namespace GeFiCa { class RhoPhi; class Segmented; }
/**
 * 2D cylindrical coordinates in rho-phi plane.
 */
class GeFiCa::RhoPhi : public Grid, public TNamed
{
   public:
      RhoPhi(size_t n1=101, size_t n2=180) : Grid(n1, n2),
      TNamed("rhophi", "2D cylindrical coordinates in rho-phi") {};

      void GetBoundaryConditionFrom(Detector &detector);
      double GetC();

      ClassDef(RhoPhi,1);
   protected:
      void OverRelaxAt(size_t idx);
      void CalculateE() ; // Fixme
};
#endif


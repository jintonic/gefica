#ifndef GeFiCa_RhoPhi
#define GeFiCa_RhoPhi

#include "Grid.h"
namespace GeFiCa { class RhoPhi; class Segmented; }
/**
 * 2D cylindrical coordinates in rho-phi plane.
 */
class GeFiCa::RhoPhi : public Grid
{
   public:
      RhoPhi(size_t n1=101, size_t n2=180) : Grid(n1, n2) {
         fName="rhophi"; fTitle="2D cylindrical coordinates in rho-phi"; }
      void SetupWith(Detector &detector);
      double GetC();
   protected:
      void OverRelaxAt(size_t idx);
      void CalculateE() ; // Fixme
      ClassDef(RhoPhi,1);
};
#endif


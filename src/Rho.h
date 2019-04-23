#ifndef GeFiCa_Rho
#define GeFiCa_Rho

#include <TNamed.h>
#include "Grid.h"
namespace GeFiCa { class Rho; }
/**
 * 1D cylindrical coordinate.
 */
class GeFiCa::Rho : public Grid, public TNamed
{
   public:
      Rho(size_t n1=101) : Grid(n1),
      TNamed("rho", "1D cylindrical coordinates") {};

      void GetBoundaryConditionFrom(Detector &detector);
      void SolveAnalytically();
      double GetC();

      ClassDef(Rho,1);
   protected:
      void OverRelaxAt(size_t idx); 
};
#endif

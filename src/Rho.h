#ifndef GeFiCa_Rho
#define GeFiCa_Rho

#include "Grid.h"
namespace GeFiCa { class Rho; }
/**
 * 1D cylindrical coordinate.
 */
class GeFiCa::Rho : public Grid
{
   public:
      Rho(size_t n1=101) :
         Grid(n1) { fName="rho"; fTitle="1D cylindrical coordinates"; }
      void SetupWith(Detector &detector);
      void SolveAnalytically();
      double GetC();
   protected:
      void OverRelaxAt(size_t idx); 
      ClassDef(Rho,1);
};
#endif

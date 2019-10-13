#ifndef GeFiCa_X
#define GeFiCa_X

#include "Grid.h"
namespace GeFiCa { class X; }
/**
 * 1D Cartesian coordinate.
 */
class GeFiCa::X : public Grid
{
   public:
      X(size_t n1=101) :
         Grid(n1) { fName="x"; fTitle="1D Cartesian coordinate"; }
      void SetupWith(Detector &detector);
      void SolveAnalytically();
      double GetC();
   protected:
      void OverRelaxAt(size_t idx); 
      ClassDef(X,1);
};
#endif

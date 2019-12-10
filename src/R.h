#ifndef GeFiCa_R
#define GeFiCa_R

#include "Grid.h"
namespace GeFiCa { class R; }
/**
 * 1D spherical coordinate.
 */
class GeFiCa::R : public Grid
{
   public:
      R(size_t n1=101) :
         Grid(n1) { fName="r"; fTitle="1D spherical coordinate"; }
      void SetupWith(Detector &detector);
      void SolveAnalytically();
      double GetC();
   protected:
      virtual void OverRelaxAt(size_t idx);
      ClassDef(R,1);
};
#endif


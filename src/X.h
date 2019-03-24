#ifndef GeFiCa_X
#define GeFiCa_X

#include "Grid.h"

namespace GeFiCa { class Grid; class X; }
/**
 * 1D Cartesian coordinate.
 */
class GeFiCa::X : public GeFiCa::Grid
{
   public:
      X(size_t nx=101) : Grid(nx) { SetName("x"); SetTitle("1D coordinate"); }
      ~X();
 
      ClassDef(X,1);

   protected:
      /**
       * Calculate fields at @param idx.
       */
      void OverRelaxAt(size_t idx); 
};
#endif


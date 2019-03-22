#ifndef GeFiCa_X
#define GeFiCa_X

#include <TNamed.h>
#include "Grid.h"

class TF3; class TTree; class TGraph;

namespace GeFiCa { class Grid; class X; }

/**
 * 1D Cartesian coordinate.
 */
class GeFiCa::X : public GeFiCa::Grid, public TNamed 
{
   public:
      X(size_t nx=101) : Grid(nx), TName("x", "1D coordinate") {};
      virtual ~X();
 
      ClassDef(X,1);

   protected:
      /**
       * Calculate fields at @param idx using SOR2.
       */
      virtual void OverRelaxAt(size_t idx); 
};
#endif


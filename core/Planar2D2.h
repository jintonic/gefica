#ifndef GeFiCa_PLANAR2D2_H
#define GeFiCa_PLANAR2D2_H

#include "XY.h"

namespace GeFiCa { class Planar2D2; }

class GeFiCa::Planar2D2 : public GeFiCa::XY
{
   public :
      double XUpperBound,XLowerBound,YUpperBound,YLowerBound; 

   public :
      Planar2D2(int ix,int iy) : XY(ix,iy),XUpperBound(1),XLowerBound(0),YUpperBound(1),YLowerBound(0){};
      void Initialize();
      bool CalculateField(EMethod method=kSOR2);

      ClassDef(Planar2D2, 1);
};

#endif

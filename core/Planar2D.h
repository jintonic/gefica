#ifndef GEFICA_PLANAR2D_H
#define GEFICA_PLANAR2D_H

#include "XY.h"

namespace GEFICA { class Planar2D; }

class GEFICA::Planar2D : public GEFICA::XY
{
   public :
      double XUpperBound,XLowerBound,YUpperBound,YLowerBound; 

      double cathode_voltage,annode_voltage;
   public :
      Planar2D(int ix,int iy) : XY(ix,iy) {};
      void initialize();
      bool CalculateField(EMethod method=kRK2);

      ClassDef(Planar2D, 1);
};

#endif

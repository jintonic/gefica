#ifndef GeFiCa_PLANAR2D_H
#define GeFiCa_PLANAR2D_H

#include "XY.h"

namespace GeFiCa { class Planar2D; }

class GeFiCa::Planar2D : public GeFiCa::XY
{
   public :
      double XUpperBound,XLowerBound,YUpperBound,YLowerBound; 

      double cathode_voltage,annode_voltage;
   public :
      Planar2D(int ix,int iy) : XY(ix,iy),XUpperBound(10),XLowerBound(0),YUpperBound(10),YLowerBound(0), cathode_voltage(2000),annode_voltage(0){};
      void Initialize();
      bool CalculateField(EMethod method=kSOR2);

      ClassDef(Planar2D, 1);
};

#endif

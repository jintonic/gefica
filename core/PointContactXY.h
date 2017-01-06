#ifndef GeFiCa_POINTCONTACTXY_H
#define GeFiCa_POINTCONTACTXY_H

#include "XY.h"

namespace GeFiCa { class PointContactXY; }

class GeFiCa::PointContactXY : public GeFiCa::XY
{
  public:
    double XUpperBound,XLowerBound,YUpperBound,YLowerBound,PointBegin,PointEnd;//bounds for X and Y and point start and end
   public :
     PointContactXY(int ix,int iy) : XY(ix,iy), XUpperBound(1),XLowerBound(0),YUpperBound(1),YLowerBound(0), PointBegin(0.4),PointEnd(0.6){};

     void Initialize();
     bool CalculateField(EMethod method=kSOR2);

     ClassDef(PointContactXY,1);
};

#endif

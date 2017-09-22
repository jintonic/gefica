#ifndef GeFiCa_POINTCONTACTRZ_H
#define GeFiCa_POINTCONTACTRZ_H

#include "RZ.h"

namespace GeFiCa { class PointContactRZ; }

class GeFiCa::PointContactRZ : public GeFiCa::RZ
{
  public:
    double XUpperBound,XLowerBound,YUpperBound,YLowerBound,PointBegin,PointEnd;//bounds for X and Y and point start and end
   public :
     PointContactRZ(int ix,int iy) : RZ(ix,iy), XUpperBound(1),XLowerBound(0),YUpperBound(1),YLowerBound(0), PointBegin(0.4),PointEnd(0.6){};

     void Initialize();
     bool CalculateField(EMethod method=kSOR2);

     ClassDef(PointContactRZ,1);
};

#endif

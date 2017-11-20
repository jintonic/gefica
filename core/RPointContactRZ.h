#ifndef GeFiCa_RPOINTCONTACTRZ_H
#define GeFiCa_RPOINTCONTACTRZ_H

#include "RZ.h"

namespace GeFiCa { class RPointContactRZ; }

class GeFiCa::RPointContactRZ : public GeFiCa::RZ
{
  public:
    double RUpperBound,RLowerBound,ZUpperBound,ZLowerBound,PointBegin,PointEnd,RHole,DHole;//bounds for X and Y and point start and end
   public :
     RPointContactRZ(int ix,int iy) : RZ(ix,iy), RUpperBound(1),RLowerBound(0),ZUpperBound(1),ZLowerBound(0), PointBegin(0.4),PointEnd(0.6),RHole(0.3),DHole(0.2){};

     void Initialize();
     bool CalculateField(EMethod method=kSOR2);

     ClassDef(RPointContactRZ,1);
};

#endif

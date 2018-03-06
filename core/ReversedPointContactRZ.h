#ifndef GeFiCa_RPOINTCONTACTRZ_H
#define GeFiCa_RPOINTCONTACTRZ_H

#include "RZ.h"

namespace GeFiCa { class ReversedPointContactRZ; }

class GeFiCa::ReversedPointContactRZ : public GeFiCa::RZ
{
  public:
    double RUpperBound, RLowerBound, ZUpperBound, ZLowerBound, PointBegin,
           PointEnd, InnerRadiusHole, OutterRadiusHole, removedConnorradius,
           removedConnorheight, DHole; // bounds for X and Y and point start and end

   public :
     ReversedPointContactRZ(int ix,int iy) : RZ(ix,iy), RUpperBound(1),
     RLowerBound(0), ZUpperBound(1), ZLowerBound(0), PointBegin(0.4),
     PointEnd(0.6), InnerRadiusHole(0.3), OutterRadiusHole(0.5),
     removedConnorradius(0.2), removedConnorheight(0.3), DHole(0.2) {};

     void Initialize();
     bool CalculateField(EMethod method=kSOR2);

     ClassDef(ReversedPointContactRZ,1);
};

#endif

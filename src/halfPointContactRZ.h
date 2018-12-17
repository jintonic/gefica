#ifndef GeFiCa_HALFPOINTCONTACTRZ_H
#define GeFiCa_HALFPOINTCONTACTRZ_H

#include "RZ.h"

namespace GeFiCa { class halfPointContactRZ; }

class GeFiCa::halfPointContactRZ : public GeFiCa::RZ
{
  public:
    double Radius,ZUpperBound,ZLowerBound,PointR,PointDepth;//bounds for X and Y and point start and end
   public :
     halfPointContactRZ(int ix=0,int iy=0) : RZ(ix,iy),
     Radius(1),ZUpperBound(1),ZLowerBound(0), PointR(0.4),PointDepth(0.2){};

     void Initialize();
     bool CalculatePotential(EMethod method=kSOR2);

     ClassDef(halfPointContactRZ,1);
};

#endif

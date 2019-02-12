#ifndef GeFiCa_RPOINTCONTACTRZ_H
#define GeFiCa_RPOINTCONTACTRZ_H

#include "RZ.h"

namespace GeFiCa { class ReversedCoaxialRZ; }

class GeFiCa::ReversedCoaxialRZ : public GeFiCa::RZ
{
  public:
    double Radius, Z,Z0, PointContactR, PointContactZ,HoleInnerR, HoleOutterR, ConnorLength,
           ConnorZ, HoleZ; // bounds for X and Y and point start and end

   public :
     ReversedCoaxialRZ(int ix,int iy) : RZ(ix,iy), Radius(1),
     Z(1), Z0(0), PointContactR(0.4),PointContactZ(0), HoleInnerR(0.3), HoleOutterR(0.5),
     ConnorLength(0.2), ConnorZ(0.3), HoleZ(0.2) {};

     void Initialize();
     bool CalculatePotential(EMethod method=kSOR2);
     void SetupBoundary();

     ClassDef(ReversedCoaxialRZ,1);
};

#endif

#ifndef GeFiCa_RPOINTCONTACTRZ_H
#define GeFiCa_RPOINTCONTACTRZ_H

#include "RZ.h"

namespace GeFiCa { class ReversedCoaxialRZ; }

class GeFiCa::ReversedCoaxialRZ : public GeFiCa::RZ
{
  public:
    double Radius, Z,Z0, Rpc, Zpc,InnerRadiusHole, OutterRadiusHole, removedConnorradius,
           removedConnorheight, DHole; // bounds for X and Y and point start and end

   public :
     ReversedCoaxialRZ(int ix,int iy) : RZ(ix,iy), Radius(1),
     Z(1), Z0(0), Rpc(0.4),Zpc(0), InnerRadiusHole(0.3), OutterRadiusHole(0.5),
     removedConnorradius(0.2), removedConnorheight(0.3), DHole(0.2) {};

     void Initialize();
     bool CalculatePotential(EMethod method=kSOR2);

     ClassDef(ReversedCoaxialRZ,1);
};

#endif

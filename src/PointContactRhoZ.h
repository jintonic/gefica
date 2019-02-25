#ifndef GeFiCa_HALFPOINTCONTACTRhoZ_H
#define GeFiCa_HALFPOINTCONTACTRhoZ_H

#include "RhoZ.h"

namespace GeFiCa { class PointContactRhoZ; }

/**
 * Grid setup for 2D point contact detectors.
 * The grid is setup in [0, Radius] and [0, Z0].
 */
class GeFiCa::PointContactRhoZ : public GeFiCa::RhoZ
{
  public:
    double Radius,Z,Z0,PointContactR,PointContactZ;//bounds for X and Y and point start and end
   public :
     PointContactRhoZ(int ix=0,int iy=0) : RhoZ(ix,iy),
     Radius(1),Z(1),Z0(0), PointContactR(0.4),PointContactZ(0.2){};

     void Initialize();
     void BoundaryOnPointcontact();

     ClassDef(PointContactRhoZ,1);
};

#endif

#ifndef GeFiCa_HALFPOINTCONTACTRZ_H
#define GeFiCa_HALFPOINTCONTACTRZ_H

#include "RZ.h"

namespace GeFiCa { class PointContactRZ; }

/**
 * Grid setup for 2D point contact detectors.
 * The grid is setup in [0, Radius] and [0, Z0].
 */
class GeFiCa::PointContactRZ : public GeFiCa::RZ
{
  public:
    double Radius,Z,Z0,PointContactR,PointContactZ;//bounds for X and Y and point start and end
   public :
     PointContactRZ(int ix=0,int iy=0) : RZ(ix,iy),
     Radius(1),Z(1),Z0(0), PointContactR(0.4),PointContactZ(0.2){};

     void Initialize();
     void BoundaryOnPointcontact();

     ClassDef(PointContactRZ,1);
};

#endif

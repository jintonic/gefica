#ifndef GeFiCa_POINTCONTACTRZ_H
#define GeFiCa_POINTCONTACTRZ_H

#include "RZ.h"

namespace GeFiCa { class PointContactRZ; }

class GeFiCa::PointContactRZ : public GeFiCa::RZ
{
   public:
      double Radius,ZUpperBound,ZLowerBound,PointR,PointDepth;//bounds for X and Y and point start and end
   public :
      PointContactRZ(int ix, int iy) : RZ(ix, iy),
      Radius(1),ZUpperBound(1),ZLowerBound(0), PointR(0.4),PointDepth(0.2){};

      /**
       * Assign initial voltage values to grid points.
       */
      void Initialize();
      /**
       * Move grids close to point contact boundary to the boundary.
       */
      void BoundaryOnPointcontact();
      bool CalculatePotential(EMethod method=kSOR2);

      ClassDef(PointContactRZ,1);

   protected:
      void SOR2(int idx, bool elec);
};

#endif

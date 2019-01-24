#ifndef GeFiCa_POINTCONTACTDZ_H
#define GeFiCa_POINTCONTACTDZ_H

#include "RZ.h"

namespace GeFiCa { class PointContactDZ; }

class GeFiCa::PointContactDZ : public GeFiCa::RZ
{
   public:
      double Radius,ZUpperBound,ZLowerBound,PointR,PointDepth,ContactInnerR;//bounds for X and Y and point start and end
 
      PointContactDZ(int ix, int iy) : RZ(ix, iy),
      Radius(1),ZUpperBound(1),ZLowerBound(0), PointR(0.4),PointDepth(0.2),ContactInnerR(1){};

      /**
       * Assign initial voltage values to grid points.
       */
      void Initialize();
      /**
       * Move grids close to point contact boundary to the boundary.
       */
      void BoundaryOnPointcontact();
      void BoundaryonWarpAround();
      bool CalculatePotential(EMethod method=kSOR2);
      bool SaveFieldasFieldgen(const char * fout);

      ClassDef(PointContactDZ,1);

   protected:
      bool CalculateField(int idx);
};

#endif

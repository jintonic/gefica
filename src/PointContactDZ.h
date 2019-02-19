#ifndef GeFiCa_POINTCONTACTDZ_H
#define GeFiCa_POINTCONTACTDZ_H

#include "RZ.h"

namespace GeFiCa { class PointContactDZ; }

/**
 * Grid setup for 2D point contact detectors.
 * The grid is setup in [-Radius, Radius] and [0, Z0].
 */
class GeFiCa::PointContactDZ : public GeFiCa::RZ
{
   public:
      double Z0, HoleInnerR, HoleOutterR, ConnorLength,
           ConnorZ, HoleZ; // bounds for X and Y and point start and end
      double Radius,Z,PointContactR,PointContactZ,WrapArroundR,TaperLength,TaperZ;//bounds for X and Y and point start and end
 
     
      PointContactDZ(int ix=101, int iy=101) : RZ(ix, iy),
      Z0(0),  HoleInnerR(0.3), HoleOutterR(0.5),
      ConnorLength(0.2), ConnorZ(0.3), HoleZ(0.2), 
      Radius(1),Z(1), PointContactR(0.4),PointContactZ(0.2),WrapArroundR(1),TaperLength(0.1),TaperZ(0.1){};

      /**
       * Assign initial voltage values to grid points.
       */
      void Initialize();
      /**
       * Move grids close to point contact boundary to the boundary.
       */
      void SetBoundary();
      bool SaveFieldAsFieldgen(const char * fout);

      ClassDef(PointContactDZ,1);

   protected:
      bool CalculateField(int idx);
};

#endif

#ifndef GeFiCa_POINTCONTACTDZ_H
#define GeFiCa_POINTCONTACTDZ_H

#include "RhoZ.h"

namespace GeFiCa { class PointContactDZ; }

/**
 * Grid setup for 2D point contact detectors.
 * The grid is setup in [-Radius, Radius] and [0, Z].
 */
class GeFiCa::PointContactDZ : public GeFiCa::RhoZ
{
   public:
      double Z0, HoleInnerR, HoleOutterR, ConnorLength,
           ConnorZ, HoleZ; // bounds for X and Y and point start and end
      double Radius,Z,PointContactR,PointContactZ,WrapArroundR,TaperLength,TaperZ;//bounds for X and Y and point start and end
 
      /**
       * Default constructor.
       */
      PointContactDZ(int ix=100, int iy=101, const char *name="pcdz",
            const char *title="2D point contact detector");

      virtual void Initialize();
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

#ifndef GeFiCa_POINTCONTACTDZ_H
#define GeFiCa_POINTCONTACTDZ_H

#include "RZ.h"

namespace GeFiCa { class PointContactDZ; }

class GeFiCa::PointContactDZ : public GeFiCa::RZ
{
   public:
      double Radius,Z,PointContactR,PointContactZ,WrapArroundR,TaperLength,TaperZ;//bounds for X and Y and point start and end
 
      PointContactDZ(int ix=101, int iy=101) : RZ(ix, iy),
      Radius(1),Z(1), PointContactR(0.4),PointContactZ(0.2),WrapArroundR(1),TaperLength(0),TaperZ(0),Z0(0){};

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
      double Z0;
      bool CalculateField(int idx);
};

#endif

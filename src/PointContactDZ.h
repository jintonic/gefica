#ifndef GeFiCa_POINTCONTACTDZ_H
#define GeFiCa_POINTCONTACTDZ_H

#include "RZ.h"

namespace GeFiCa { class PointContactDZ; }

class GeFiCa::PointContactDZ : public GeFiCa::RZ
{
   public:
      double Radius,Z,Rpc,Zpc,RwrapArround,TaperLength;//bounds for X and Y and point start and end
 
      PointContactDZ(int ix, int iy) : RZ(ix, iy),
      Radius(1),Z(1), Rpc(0.4),Zpc(0.2),RwrapArround(1),TaperLength(0),Z0(0){};

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
      bool SaveFieldAsFieldgen(const char * fout);

      ClassDef(PointContactDZ,1);

   protected:
      double Z0;
      bool CalculateField(int idx);
};

#endif

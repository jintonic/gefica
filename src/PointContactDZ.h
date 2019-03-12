#ifndef GeFiCa_POINTCONTACTDZ_H
#define GeFiCa_POINTCONTACTDZ_H

#include "RhoZ.h"

namespace GeFiCa { class PointContactDZ; }

/**
 * Grid setup for 2D point contact detectors.
 * The grid is set in [-Radius, Radius] and [0, Height].
 */
class GeFiCa::PointContactDZ : public GeFiCa::RhoZ
{
   public:
      double Height; ///< Height of crystal
      double Radius; ///< Radius of crystal

      double PointContactH; ///< Height of point contact
      double PointContactR; ///< Radius of point contact

      double HoleH; ///< Depth of bore hole
      double HoleR; ///< radius of bore hole
      double HoleTaperW; ///< width of bore hole taper
      double HoleTaperH; ///< height of bore hole taper

      double TaperW; ///< Width of taper (point contact side)
      double TaperH; ///<Height of taper (point contact side)

      double CornerW; ///< Width of taper (bore side)
      double CornerH; ///<Height of taper (bore side)

      double WrapArroundR; ///< Inner radius of outer contact 
      double GrooveW; ///< Width of Groove 
      double GrooveH; ///< Height of Groove 

      /**
       * Default constructor.
       */
      PointContactDZ(int nd=300, ///< [in] number of points across diameter
            int nz=301, ///< [in] number of points along height
            const char *name="pcdz", ///< [in] name of the class object created
            const char *title="2D point contact detector");

      void Export2fieldgen(const char *output="ppc.dat");

      ClassDef(PointContactDZ,1);

   protected:
      virtual void InitializeGrid();
      virtual void SetGridImpurity();
      /**
       * Move grids close to point contact boundary to the boundary.
       */
      void SetBoundary();
};

#endif

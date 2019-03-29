#ifndef GeFiCa_PointContact
#define GeFiCa_PointContact

#include "Detector.h"
namespace GeFiCa { class PointContact; }
/**
 * Configuration of point contact detectors.
 */
class GeFiCa::PointContact : public Detector
{
   public:
      double Radius; ///< Radius of crystal

      double PointContactR; ///< Radius of point contact
      double PointContactH; ///< Height of point contact

      double BoreH; ///< Depth of bore hole
      double BoreR; ///< radius of bore hole
      double BoreTaperW; ///< width of bore hole taper
      double BoreTaperH; ///< height of bore hole taper

      double TaperW; ///< Width of taper (point contact side)
      double TaperH; ///<Height of taper (point contact side)

      double CornerW; ///< Width of taper (bore side)
      double CornerH; ///<Height of taper (bore side)

      double WrapAroundR; ///< Inner radius of outer contact 
      double GrooveW; ///< Width of Groove 
      double GrooveH; ///< Height of Groove 

      PointContact(const char *name="pc",
            const char *title="point-contact detector");
      void CheckConfigurations();
      ClassDef(PointContact,1);
};
#endif


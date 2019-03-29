#ifndef GeFiCa_Hemispherical
#define GeFiCa_Hemispherical

#include "Detector.h"
namespace GeFiCa { class Hemispherical; }
/**
 * Configuration of hemispherical detectors.
 */
class GeFiCa::Hemispherical : public Detector
{
  public:
      double PointContactR; ///< radius of point contact
      double PointContactH; ///< height of point contact

      Hemispherical(const char *name="hs",
            const char *title="hemispherical detector");
      void CheckConfigurations();
      ClassDef(Hemispherical, 1);
};
#endif


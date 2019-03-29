#ifndef GeFiCa_TrueCoaxial
#define GeFiCa_TrueCoaxial

#include "Detector.h"
namespace GeFiCa { class TrueCoaxial; }
/**
 * Configuration of true coaxial detectors.
 */
class GeFiCa::TrueCoaxial : public Detector
{
   public :
      double Radius; ///< radius of the detector
      double BoreR; ///< radius of the bore

      TrueCoaxial(const char *name="tc",
            const char *title="true coaxial detector");
      void CheckConfigurations();
      ClassDef(TrueCoaxial, 1);
};
#endif

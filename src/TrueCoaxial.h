#ifndef GeFiCa_TrueCoaxial
#define GeFiCa_TrueCoaxial

#include "Detector.h"
namespace GeFiCa { class TrueCoaxial; }
/**
 * Configuration of a true coaxial detector.
 */
class GeFiCa::TrueCoaxial : public Detector
{
   public :
      double Radius; ///< radius of the detector
      double BoreR; ///< radius of the bore

      TrueCoaxial(const char *name="tc",
            const char *title="true coaxial detector");
      ClassDef(TrueCoaxial, 1);
};
#endif

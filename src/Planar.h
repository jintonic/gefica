#ifndef GeFiCa_Planar
#define GeFiCa_Planar

#include "Detector.h"
namespace GeFiCa { class Planar; }
/**
 * Configuration of planar detectors.
 */
class GeFiCa::Planar : public Detector
{
   public :
      double Width; ///< width of a planar detector
      double Depth; ///< depth of a planar detector

      Planar(const char *name="planar", const char *title="planar detector");
      void CheckConfigurations();
      ClassDef(Planar, 1);
};
#endif

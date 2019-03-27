#ifndef GeFiCa_Planar
#define GeFiCa_Planar

#include "Detector.h"
namespace GeFiCa { class Planar; }
/**
 * Configuration of a planar detector.
 */
class GeFiCa::Planar : public GeFiCa::Detector
{
   public :
      double Width; ///< width of a planar detector
      double Depth; ///< depth of a planar detector

      Planar(const char *name="planar", const char *title="planar detector");
      ClassDef(Planar, 1);
};
#endif

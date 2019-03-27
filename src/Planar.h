#ifndef GeFiCa_Planar
#define GeFiCa_Planar

#include <TNamed.h>
#include "Detector.h"
namespace GeFiCa { class Planar; }
/**
 * Configuration of a planar detector.
 */
class GeFiCa::Planar : public GeFiCa::Detector, public TNamed
{
   public :
      double Width; ///< width of a planar detector
      double Depth; ///< depth of a planar detector

      Planar();
      ClassDef(Planar, 1);
};
#endif

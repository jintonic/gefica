#ifndef GeFiCa_Planar
#define GeFiCa_Planar

#include <TObject.h>
#include "Detector.h"
namespace GeFiCa { class Planar; }
/**
 * Configuration of a planar detector.
 */
class GeFiCa::Planar : public GeFiCa::Detector, public TObject
{
   public :
      double Width; ///< width of a planar detector

      Planar();
      ClassDef(Planar, 1);
};
#endif

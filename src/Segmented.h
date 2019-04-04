#ifndef GeFiCa_Segmented
#define GeFiCa_Segmented

#include "Detector.h"

namespace GeFiCa { class Segmented; }

/**
 * Configuration of segmented true coaxial detectors.
 */
class GeFiCa::Segmented: public Detector
{
   public:
      double Radius; ///< radius of the crystal
      double BoreR; ///< radius of the bore hole
      size_t Nphi; ///< total number of segments in phi
      size_t Nz; ///< total number of segments in z
      size_t SegmentId; ///< segment Id in [0, Nphi*Nz]
 
      Segmented(const char *name="sip",
            const char *title="segmented detector");
      void CheckConfigurations();
      void Draw(Option_t* option="");
      ClassDef(Segmented,1);
};
#endif


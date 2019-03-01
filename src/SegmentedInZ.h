#ifndef GeFiCa_SEGMENTEDINZ_H
#define GeFiCa_SEGMENTEDINZ_H

#include "RhoZ.h"

namespace GeFiCa { class SegmentedInZ; }

/**
 * Grid setup for 2D segmented true coaxial detectors.
 */
class GeFiCa::SegmentedInZ : public GeFiCa::RhoZ
{
   public :
      double InnerR;///< inner radius
      double OuterR;///< outer radius
      double Z; ///< height of detector along z
      int Nseg; ///< number of segments 
      int SegmentId; ///< segment Id

      /**
       * Default constructor.
       */
      SegmentedInZ(int nr=300, ///< [in] number of grid points along r
            int nt=300, ///< [in] number of grid points along theta
            const char *name="siz", ///< [in] name of object
            const char *title="2D true coaxial detector segmented in z");

      ClassDef(SegmentedInZ,1);

   protected:
      virtual void Initialize(); 
};
#endif


#ifndef GeFiCa_SIEGFRIEDInPhi_H
#define GeFiCa_SIEGFRIEDInPhi_H

#include "RhoPhi.h"
#include <TMath.h> 

namespace GeFiCa { class SegmentedInPhi; }

/**
 * Grid setup for 2D true coaxial detectors segmented in phi.
 */
class GeFiCa::SegmentedInPhi : public GeFiCa::RhoPhi
{
   public:
      double OuterR; ///< outer radius
      double InnerR; ///< inner radius
      int Nseg; ///< total number of segments
      int SegmentId; ///< segment Id
 
      /**
       * Default constructor
       */
      SegmentedInPhi(int nr=360, ///< [in] number of grid points along rho
            int np=301, ///< [in] number of grid points along phi
            const char *name="sip", ///< [in] name of object
            const char *title="2D true coaxial detector segmented in phi");

      virtual void Initialize();

      ClassDef(SegmentedInPhi,1);
};

#endif

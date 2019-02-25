#ifndef GeFiCa_SIEGFRIEDInPhi_H
#define GeFiCa_SIEGFRIEDInPhi_H

#include "RhoPhi.h"
#include <TMath.h> 

namespace GeFiCa { class SegmentedInPhi; }

/**
 * Grid setup for InPhi segmented true coaxial detectors.
 */
class GeFiCa::SegmentedInPhi : public GeFiCa::RhoPhi
{
   public:
      double RUpperBound,RLowerBound;
      int SegmentNum;//total segment number
      int SegmentID;
 
      SegmentedInPhi(int ix=360, int iy=301) : RhoPhi(ix, iy),
         RUpperBound(2.5),RLowerBound(0.5),SegmentNum(6) ,SegmentID(1) {};

      virtual void Initialize();

      ClassDef(SegmentedInPhi,1);
};

#endif

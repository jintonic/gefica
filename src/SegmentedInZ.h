#ifndef GeFiCa_SEGMENTEDINZ_H
#define GeFiCa_SEGMENTEDINZ_H

#include "RhoZ.h"

namespace GeFiCa { class SegmentedInZ; }

/**
 * Grid setup for 2D true coaxial detectors.
 */
class GeFiCa::SegmentedInZ : public GeFiCa::RhoZ
{
   public :
      double InnerR,OuterR,Height;
      int SegmentNum;//total segment number
      int SegmentID;

      SegmentedInZ(int r=300,int O=300) : RhoZ(r, O) ,InnerRadius(0.5),OuterRadius(3),Height(3),SegmentNum(6) ,SegmentID(1) {};

      virtual void Initialize(); 

      ClassDef(SegmentedInZ,1);
};
#endif


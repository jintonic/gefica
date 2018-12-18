#ifndef GeFiCa_SIEGFRIED2D_H
#define GeFiCa_SIEGFRIED2D_H

#include "RhoPhi.h"
#include <TMath.h> 

namespace GeFiCa { class Segmented2D; }

class GeFiCa::Segmented2D : public GeFiCa::RhoPhi
{
   public:
      double RUpperBound,RLowerBound;
      int SegmentNum;//total segment number
 
      Segmented2D(int ix=360, int iy=301) : RhoPhi(ix, iy),
         RUpperBound(2.5),RLowerBound(0.5),SegmentNum(6) {};

      /**
       * Assign initial voltage values to grid points.
       */
      void Initialize(int SegmentID);
      bool CalculatePotential(EMethod method=kSOR2,int SegmentID=1);

      ClassDef(Segmented2D,1);
};

#endif

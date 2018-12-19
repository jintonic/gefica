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
      int SegmentID;
 
      Segmented2D(int ix=360, int iy=301) : RhoPhi(ix, iy),
         RUpperBound(2.5),RLowerBound(0.5),SegmentNum(6) ,SegmentID(1) {};

      /**
       * Assign initial voltage values to grid points.
       */
      void Initialize();
      bool CalculatePotential(EMethod method=kSOR2);

      ClassDef(Segmented2D,1);
};

#endif

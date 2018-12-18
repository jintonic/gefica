#ifndef GeFiCa_SIEGFRIED2D_H
#define GeFiCa_SIEGFRIED2D_H

#include "RhoPhi.h"
#include <TMath.h> 

namespace GeFiCa { class Siegfried2D; }

class GeFiCa::Siegfried2D : public GeFiCa::RhoPhi
{
   public:
      double RUpperBound,RLowerBound,SegmentSize;//bounds for X and Y and point start and end
 
      Siegfried2D(int ix, int iy) : RhoPhi(ix, iy),
      RUpperBound(2.5),RLowerBound(0.5),SegmentSize(TMath::Pi()/3) {};

      /**
       * Assign initial voltage values to grid points.
       */
      void Initialize();
      bool CalculatePotential(EMethod method=kSOR2);

      ClassDef(Siegfried2D,1);
};

#endif

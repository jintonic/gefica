#ifndef GeFiCa_RHOPHI_H
#define GeFiCa_RHOPHI_H

#include "XY.h"

namespace GeFiCa { class RhoPhi; }

/**
 * 2D cylindrical coordinates.
 */
class GeFiCa::RhoPhi : public GeFiCa::XY
{
   public:
      RhoPhi(int nRho=101, int nPhi=101): XY(nRho,nPhi) {}; 

      ClassDef(RhoPhi,1);

   protected:
      void DoSOR2(int idx);
      virtual double GetData(double tarx,double tary,EOutput output);
};
#endif


#ifndef GeFiCa_RZ_H
#define GeFiCa_RZ_H

#include "XY.h"
class TF2;

namespace GeFiCa { class RZ; }

/**
 *2D Field based on X
 * it will create a 2D Field using cartesian coordinate and calculate the result
 */
class GeFiCa::RZ : public GeFiCa::XY
{
   public:
      /**
       * RZ is a constructor with default values n1,n2 = 101
       */
      RZ(unsigned short n1=101, unsigned short n2=101):XY(n1,n2){};

      double GetCapacitance();

      ClassDef(RZ,1);

   protected:
      virtual void SOR2(int idx,bool elec);
};
#endif 

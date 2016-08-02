#ifndef GEFICA_RHOPHI_H
#define GEFICA_RHOPHI_H

#include "XY.h"

namespace GEFICA { class RhoPhi; }

class GEFICA::RhoPhi : public GEFICA::XY
{
   public:
      RhoPhi(unsigned short r=101, unsigned short phi=101): XY(r,phi) {} 
      virtual ~RhoPhi() {};

      virtual void CreateGridWithFixedStepLength(double steplength); 

      virtual double GetData(double tarx,double tary,int thing);
      void RK2(int idx,bool elec);

      ClassDef(RhoPhi,1);
};

#endif

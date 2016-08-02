#ifndef GEFICA_POLAR_H
#define GEFICA_POLAR_H

#include "XY.h"

namespace GEFICA { class Polar; }

class GEFICA::Polar : public GEFICA::XY
{
   public:
      Polar(unsigned short r=101, unsigned short phi=101): XY(r,phi) {} 
      virtual ~Polar() {};

      virtual void CreateGridWithFixedStepLength(double steplength); 

      virtual double GetData(double tarx,double tary,int thing);
      void RK2(int idx,bool elec);

      ClassDef(Polar,1);
};

#endif

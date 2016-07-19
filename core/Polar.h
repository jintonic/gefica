#ifndef GEFICA_POLAR_H
#define GEFICA_POLAR_H

#include "XY.h"

namespace GEFICA { class Polar; }

class GEFICA::Polar : public GEFICA::XY
{
   public:
      Polar(unsigned short x=0, unsigned short y=0): XY(x,y) {} 
      virtual ~Polar();

      virtual void Create(double steplength); 

      virtual double GetData(double tarx,double tary,int thing);

      ClassDef(Polar,1);
};

#endif

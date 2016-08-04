#ifndef GEFICA_SPHERIC_H
#define GEFICA_SPHERIC_H

#include "RThetaPhi.h"

namespace GEFICA { class Spheric;}

class GEFICA::Spheric: public GEFICA::RThetaPhi
{
   public :
      double r1,r0;

      Spheric(int r,int O,int a) : RThetaPhi(r, O,a ) {};
      void SetVoltage(double v0,double v1); 
      void Create(double r0,double r1);

      ClassDef(Spheric,1);
};

#endif

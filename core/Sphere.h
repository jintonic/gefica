#ifndef GEFICA_SPHERE_H
#define GEFICA_SPHERE_H

#include "RThetaPhi.h"

namespace GEFICA { class Sphere;}

class GEFICA::Sphere: public GEFICA::RThetaPhi
{
   public :
      double r1,r0;

      Sphere(int r,int O,int a) : RThetaPhi(r, O,a ) {};

      void SetVoltage(double v1,double v2);

      ClassDef(Sphere,1);
};

#endif

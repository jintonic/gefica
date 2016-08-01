#ifndef GEFICA_BALL_H
#define GEFICA_BALL_H

#include "Spherical.h"

namespace GEFICA { class ball; }

class GEFICA::ball : public GEFICA::Spherical
{
   public :
      double r1,r0;

      ball(int r,int O,int a) : Spherical(r, O,a ) {};
      void SetVoltage(double v0,double v1); 
      void Create(double r0,double r1);

      ClassDef(ball,1);
};

#endif

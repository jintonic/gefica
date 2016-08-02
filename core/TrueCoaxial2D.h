#ifndef GEFICA_TRUECOAXIAL2D_H
#define GEFICA_TRUECOAXIAL2D_H

#include "RhoPhi.h"

namespace GEFICA { class TrueCoaxial2D; }

class GEFICA::TrueCoaxial2D : public GEFICA::RhoPhi
{
   public :
      double r1,r0;

      TrueCoaxial2D(int r,int O) : RhoPhi(r, O) {};
      void SetVoltage(double v0,double v1); 
      void Create(double r0,double r1);

      ClassDef(TrueCoaxial2D,1);
};

#endif

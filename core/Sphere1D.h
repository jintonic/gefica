#ifndef GEFICA_SPHERE1D_H
#define GEFICA_SPHERE1D_H

#include "R.h"

namespace GEFICA { class Sphere1D; }

class GEFICA::Sphere1D : public GEFICA::R
{


   public :
      Sphere1D(int nx=101) : R(nx) {};
      ClassDef(Sphere1D, 1);

   protected:
      bool Analyic();
};

#endif

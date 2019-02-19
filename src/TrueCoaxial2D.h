#ifndef GeFiCa_TRUECOAXIAL2D_H
#define GeFiCa_TRUECOAXIAL2D_H

#include "RhoPhi.h"

namespace GeFiCa { class TrueCoaxial2D; }

/**
 * Grid setup for 2D true coaxial detectors.
 */
class GeFiCa::TrueCoaxial2D : public GeFiCa::RhoPhi
{
   public :
      double InnerRadius,OuterRadius;

      TrueCoaxial2D(int r,int O) : RhoPhi(r, O) ,InnerRadius(0.5),OuterRadius(3){};

      virtual void Initialize(); 

      ClassDef(TrueCoaxial2D,1);
};
#endif


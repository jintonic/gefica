#ifndef GeFiCa_TRUECOAXIAL2D2_H
#define GeFiCa_TRUECOAXIAL2D2_H

#include "RZ.h"

namespace GeFiCa { class TrueCoaxial2D2; }

class GeFiCa::TrueCoaxial2D2 : public GeFiCa::RZ
{
   public :
      double InnerRadius,OuterRadius,height;

   public:

      TrueCoaxial2D2(int r,int O) : RZ(r, O) ,InnerRadius(0.5),OuterRadius(3),height(5){};
      void Initialize(); 
      bool CalculatePotential(EMethod method=kSOR2);


      ClassDef(TrueCoaxial2D2,1);
};
#endif


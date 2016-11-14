#ifndef GeFiCa_TRUECOAXIAL2D_H
#define GeFiCa_TRUECOAXIAL2D_H

#include "RhoPhi.h"

namespace GeFiCa { class TrueCoaxial2D; }

class GeFiCa::TrueCoaxial2D : public GeFiCa::RhoPhi
{
   public :
      double InnerRadius,OuterRadius;
      double Vpos,Vneg;

   public:

      TrueCoaxial2D(int r,int O) : RhoPhi(r, O) ,InnerRadius(0.3),OuterRadius(3),Vpos(2000),Vneg(0){};
      void Initialize(); 
      bool CalculateField(EMethod method=kSOR2);


      ClassDef(TrueCoaxial2D,1);
};
#endif


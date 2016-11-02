#ifndef GeFiCa_TRUECOAXIAL2D_H
#define GeFiCa_TRUECOAXIAL2D_H

#include "RhoPhi.h"

namespace GeFiCa { class TrueCoaxial2D; }

class GeFiCa::TrueCoaxial2D : public GeFiCa::RhoPhi
{
   public :
      double InnerRadius,OuterRadius;
      double cathode_voltage,anode_voltage;

   public:

      TrueCoaxial2D(int r,int O) : RhoPhi(r, O) ,InnerRadius(0.3),OuterRadius(3),cathode_voltage(2000),anode_voltage(0){};
      void Initialize(); 
      bool CalculateField(EMethod method=kSOR2);


      ClassDef(TrueCoaxial2D,1);
};

#endif

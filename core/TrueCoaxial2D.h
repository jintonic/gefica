#ifndef GEFICA_TRUECOAXIAL2D_H
#define GEFICA_TRUECOAXIAL2D_H

#include "RhoPhi.h"

namespace GEFICA { class TrueCoaxial2D; }

class GEFICA::TrueCoaxial2D : public GEFICA::RhoPhi
{
   public :
      double InnerRadius,OuterRadius;
      double cathode_voltage,anode_voltage;

   public:

      TrueCoaxial2D(int r,int O) : RhoPhi(r, O) ,InnerRadius(0.3),OuterRadius(3),cathode_voltage(2000),anode_voltage(0){};
      void initialize(); 
      bool CalculateField(EMethod method=kRK2);


      ClassDef(TrueCoaxial2D,1);
};

#endif

#ifndef GEFICA_SPHERE1D_H
#define GEFICA_SPHERE1D_H

#include "R.h"

namespace GEFICA { class Sphere1D; }

class GEFICA::Sphere1D : public GEFICA::R
{
  public:
      double innerR,outterR; // boundary of the planar detector
      double cathode_voltage,annode_voltage;

   public :
      Sphere1D(int nx=101) : R(nx),innerR(0.3),outterR(3),cathode_voltage(2000),annode_voltage(0) {};
      ClassDef(Sphere1D, 1);
      void initialize();      
      bool CalculateField(EMethod method=kSOR2);
   protected:
      bool Analyic();
};

#endif

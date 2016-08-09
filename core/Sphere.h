#ifndef GEFICA_SPHERE_H
#define GEFICA_SPHERE_H

#include "RThetaPhi.h"

namespace GEFICA { class Sphere;}

class GEFICA::Sphere: public GEFICA::RThetaPhi
{
   public :
     double UpperBound,LowerBound; // boundary of the planar detector   
     double cathode_voltage,anode_voltage;
   public:
      Sphere(int r,int O,int a) : RThetaPhi(r, O,a ) {};

      void initialize();
      bool CalculateField(EMethod method=kRK2);

      ClassDef(Sphere,1);
};

#endif

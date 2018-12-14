#ifndef GeFiCa_SPHERE_H
#define GeFiCa_SPHERE_H

#include "RThetaPhi.h"

namespace GeFiCa { class Sphere;}

/**
 * Sphere detector in RThetaPhi coordinates.
 */
class GeFiCa::Sphere: public GeFiCa::RThetaPhi
{
   public :
      double OuterRadius,InnerRadius;
   public:
      Sphere(int r,int O,int a) : RThetaPhi(r, O, a),
      OuterRadius(3),InnerRadius(0.3) {};

      void Initialize();
      bool CalculatePotential(EMethod method=kSOR2);

      ClassDef(Sphere,1);
};
#endif


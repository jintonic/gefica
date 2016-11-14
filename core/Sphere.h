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
      double UpperBound,LowerBound;
      double Vpos,Vneg;
   public:
      Sphere(int r,int O,int a) : RThetaPhi(r, O, a),
      UpperBound(10),LowerBound(0), Vpos(2000),Vneg(0){};

      void Initialize();
      bool CalculateField(EMethod method=kSOR2);

      ClassDef(Sphere,1);
};
#endif


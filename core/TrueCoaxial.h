#ifndef TRUECOAXIAL_H
#define TRUECOAXIAL_H

#include "Detector.h"

namespace GEFICA { class TrueCoaxial; }

class GEFICA::TrueCoaxial : public GEFICA::Detector
{
   public:
      TrueCoaxial(double ri, double ro): Detector(), ri(ri), ro(ro) {};
      virtual ~TrueCoaxial() {};

      double ri, ro;
};

#endif


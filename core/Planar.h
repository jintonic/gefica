#ifndef PLANAR_H
#define PLANAR_H

#include "Detector.h"

namespace GEFICA { class Planar; }

class GEFICA::Planar : public GEFICA::Detector
{
   public:
      Planar(double h=0): Detector(), height(h) {};
      virtual ~Planar() {};

      double height;
};

#endif


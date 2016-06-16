#ifndef FIELD2D_H
#define FIELD2D_H

#include "Field.h"

namespace GEFICA { class Field2D; }

class GEFICA::Field2D : public GEFICA::Field
{
   public:
      Field2D(unsigned short n1=0, unsigned short n2=0): Field(n1), n2(n2) {};
      virtual ~Field2D() {};

      unsigned short n2; // number of steps along the 2nd axis
      std::vector<double> e2;
};

#endif


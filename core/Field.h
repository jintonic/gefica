#ifndef FIELD_H
#define FIELD_H

#include <vector>

namespace GEFICA { class Field; }

class GEFICA::Field
{
   public:
      Field(unsigned short n1=0): n1(n1) {};
      virtual ~Field() {};

      unsigned short n1; // number of steps along the 1st axis

      std::vector<double> p;
      std::vector<double> e1;
};

#endif


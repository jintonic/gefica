#ifndef DETECTOR_H
#define DETECTOR_H

class TF1;
namespace GEFICA {
   class Detector;
   class Field;
}

class GEFICA::Detector
{
   public:
      Detector(): voltage(0), impurity(0), field(0) {};
      virtual ~Detector() {};

      void UpdateField() {};

      double voltage;

      TF1 *impurity;
      Field *field;
};

#endif


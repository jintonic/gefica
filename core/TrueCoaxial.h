#ifndef FIELD2DC_H
#define TRUECOAXIAL_H

#include "Field2D.h"

namespace GEFICA { class TrueCoaxial; }

class GEFICA::TrueCoaxial : public GEFICA::Field2D
{
  public:
   TrueCoaxial(unsigned short x=0, unsigned short y=0): Field2D(x,y) {} 
    virtual ~TrueCoaxial() {delete[] E1;delete[] E2;delete[] C1; delete [] C2; delete[] P; delete[] StepNext;delete[] StepBefore;delete[] isbegin; delete[] Impurity;};


    virtual void Create(double steplength); 

    virtual double GetData(double tarx,double tary,int thing);
    ClassDef(Field2D,1);

};
#endif


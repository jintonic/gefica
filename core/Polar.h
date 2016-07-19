#ifndef POLAR_H
#define POLAR_H

#include "XY.h"

namespace GEFICA { class Polar; }

class GEFICA::Polar : public GEFICA::XY
{
  public:
   Polar(unsigned short x=0, unsigned short y=0): XY(x,y) {} 
    virtual ~Polar() {delete[] E1;delete[] E2;delete[] C1; delete [] C2; delete[] P; delete[] StepNext;delete[] StepBefore;delete[] isbegin; delete[] Impurity;};


    virtual void Create(double steplength); 

    virtual double GetData(double tarx,double tary,int thing);
    ClassDef(Polar,1);

};
#endif


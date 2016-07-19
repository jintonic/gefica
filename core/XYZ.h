#ifndef FIELD3D_H
#define FIELD3D_H

#include "Field2D.h"
#include <TF3.h>

namespace GEFICA { class Field3D; }

class GEFICA::Field3D : public GEFICA::Field2D
{
   public:
    Field3D(unsigned short x=0, unsigned short y=0,unsigned short z=0): Field2D(x,y), z(z) {n=x*y*z;}
    virtual ~Field3D() {delete[] E1;delete[] E2;delete []E3;delete[] C1; delete [] C2;delete[] C3; delete[] P; delete[] StepNext;delete[] StepBefore;delete []StepLeft;delete[] StepRight; delete []StepUp; delete [] StepDown; delete[] isbegin; delete[] Impurity;};

    unsigned short z; // number of steps along the 3nd axis
    double *E3,*C3,*StepUp,*StepDown;//left and right are for y axis
      
    virtual void Create(double steplength);
    virtual void Update(int idx); 

    virtual void Save(const char *fout=NULL);
    virtual void Load(const char *fin=NULL);
    
    virtual int FindIdx(double tarx,double tary,
	double tarz,int zbegin,int zend);
    virtual double GetData(double tarx,double tary,double tarz,int thing);
    virtual void SetImpurity(TF3 * Im);
    ClassDef(Field3D,1);
};

#endif


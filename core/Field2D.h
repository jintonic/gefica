#ifndef FIELD2D_H
#define FIELD2D_H

#include "Field.h"
#include <TF2.h>

namespace GEFICA { class Field2D; }

class GEFICA::Field2D : public GEFICA::Field
{
   public:
    Field2D(unsigned short x=0, unsigned short y=0): Field(x), y(y) {n=x*y;}
    virtual ~Field2D() {delete[] E1;delete[] E2;delete[] C1; delete [] C2; delete[] P; delete[] StepNext;delete[] StepBefore;delete[] isbegin; delete[] Impurity;};

    unsigned short y; // number of steps along the 2nd axis
    double *E2,*C2,*StepLeft,*StepRight;//left and right are for y axis
      
    virtual void Create(double steplength);
    virtual void Update(int idx); 

    virtual void Save(const char *fout=NULL);
    virtual void Load(const char *fin=NULL);
    
    virtual int FindIdx(double tarx,double tary
	,int ybegin,int yend);
    virtual double GetData(double tarx,double tary,int thing);
    virtual void SetImpurity(TF2 * Im);
    ClassDef(Field2D,1);
};

#endif


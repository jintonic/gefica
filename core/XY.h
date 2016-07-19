#ifndef GEFICA_XY_H
#define GEFICA_XY_H

#include "X.h"
#include <TF2.h>

namespace GEFICA { class XY; }

class GEFICA::XY : public GEFICA::X
{
   public:
      unsigned short y; // number of steps along the 2nd axis
      double *E2,*C2,*StepLeft,*StepRight;//left and right are for y axis

   public:
      XY(unsigned short x=0, unsigned short y=0): X(x), y(y) {n=x*y;}
      virtual ~XY();

      virtual void Create(double steplength);
      virtual void Update(int idx); 

      virtual void Save(const char *fout=NULL);
      virtual void Load(const char *fin=NULL);

      virtual double GetData(double tarx,double tary,int thing);
      virtual void SetImpurity(TF2 * Im);

      ClassDef(XY,1);

   protected:
      virtual int FindIdx(double tarx,double tary
            ,int ybegin,int yend);
};

#endif

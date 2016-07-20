#ifndef GEFICA_XY_H
#define GEFICA_XY_H

#include "X.h"
#include <TF2.h>

namespace GEFICA { class XY; }

class GEFICA::XY : public GEFICA::X
{
   public:
      unsigned short n2; // number of steps along the 2nd axis

   public:
      XY(unsigned short x=0, unsigned short y=0): X(n1), y(n2) {n=n1*n2;}
      virtual ~XY();

      virtual void Create(double steplength);
      virtual void Update(int idx); 

      virtual void Save(const char *fout=NULL);
      virtual void Load(const char *fin=NULL);

      virtual double GetData(double tarx,double tary,int thing);
      virtual void SetImpurity(TF2 * Im);

      ClassDef(XY,1);

   protected:

      double *fE2,*fC2,*fDistanceToLeft,*fDistanceToRight;//left and right are for y axis
      virtual int FindIdx(double tarx,double tary
            ,int ybegin,int yend);
};

#endif

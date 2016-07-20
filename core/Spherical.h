#ifndef GEFICA_HALFBALL_H
#define GEFICA_HALFBALL_H

#include "XYZ.h"

namespace GEFICA { class Spherical; }

class GEFICA::Spherical : public GEFICA::XYZ
{
   public:

   public:
      Spherical(unsigned short x=0, unsigned short y=0,unsigned short z=0): XYZ(x,y,z) {}
      virtual ~Spherical();

      virtual void Create(double steplength);
      virtual void Update(int idx); 

      virtual void Save(const char *fout=NULL);
      virtual void Load(const char *fin=NULL);

      virtual double GetData(double tarx,double tary,double tarz,int thing);
      virtual void SetImpurity(TF3 * Im);

      ClassDef(Spherical,1);
      
   protected:
      double *fE3,*fC3,*fDistanceToUp,*fDistanceToDown;//left and right are for y axis
      virtual int FindIdx(double tarx,double tary,
            double tarz,int zbegin,int zend);
};

#endif

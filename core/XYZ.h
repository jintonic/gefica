#ifndef GEFICA_XYZ_H
#define GEFICA_XYZ_H

class TF3;

#include "XY.h"

namespace GEFICA { class XYZ; }

class GEFICA::XYZ : public GEFICA::XY
{
   public:
      unsigned short n3; // number of steps along the 3nd axis

   public:
      XYZ(unsigned short n1=101, unsigned short n2=11,unsigned short n3=11);
      virtual ~XYZ();

      virtual void CreateGridWithFixedStepLength(double steplength);
      virtual void Update(int idx); 

      virtual void Save(const char *fout=NULL);
      virtual void Load(const char *fin=NULL);

      virtual double GetData(double tarx,double tary,double tarz,int thing);
      virtual void SetImpurity(TF3 * Im);
      
      ClassDef(XYZ,1);

   protected:

      double *fE3,*fC3,*fDistanceToUp,*fDistanceToDown;//left and right are for y axis
      virtual int FindIdx(double tarx,double tary,
            double tarz,int begin,int end);
};

#endif

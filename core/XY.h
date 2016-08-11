///////////////////////
//                   //
//XY                 //
//                   //
//2D Field based on X//
//it will create a 2D Field using cartesian coordinate and calculate the result
///////////////////////
#ifndef GEFICA_XY_H
#define GEFICA_XY_H

#include "X.h"
class TF2;

namespace GEFICA { class XY; }

class GEFICA::XY : public GEFICA::X
{
   public:
      unsigned short n2; // number of steps along the 2nd axis

   public:
      XY(unsigned short n1=101, unsigned short n2=101);
      
      virtual ~XY();

      virtual void SetStepLength(double steplength1,double steplength2);
      virtual void SOR2(int idx,bool elec); 

      virtual void SaveField(const char *fout=NULL);
      virtual void LoadField(const char *fin=NULL);

      virtual double GetData(double tarx,double tary,int thing);
      virtual double GetPotential(double tarx,double tary){return GetData(tarx,tary,1);};
      virtual double GetE1(double tarx,double tary){return GetData(tarx,tary,2);};
      virtual double GetE2(double tarx,double tary){return GetData(tarx,tary,3);};
      virtual double GetImpurity(double tarx,double tary){return GetData(tarx,tary,0);};

   
      virtual void SetImpurity(TF2 * Im);
      virtual void SetImpurity(double density);

      ClassDef(XY,1);

   protected:

      double *fE2,*fC2,*fDistanceToLeft,*fDistanceToRight;//left and right are for y axis
      virtual int FindIdx(double tarx,double tary
            ,int ybegin,int yend);
};

#endif

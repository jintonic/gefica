///////////////////////
//                   //
//XY                 //
//                   //
//2D Field based on X//
//it will create a 2D Field using cartesian coordinate and calculate the result
///////////////////////
#ifndef GeFiCa_XY_H
#define GeFiCa_XY_H

#include "X.h"
class TF2;

namespace GeFiCa { class XY; }

class GeFiCa::XY : public GeFiCa::X
{
   public:
      unsigned short n2; // number of steps along the 2nd axis

   public:
      XY(unsigned short n1=101, unsigned short n2=101);
      
      virtual ~XY();


      virtual void SaveField(const char *fout=NULL);
      virtual void LoadField(const char *fin=NULL);

      virtual double GetPotential(double x,double y){return GetData(x,y,1);};
      virtual double GetE1(double x,double y){return GetData(x,y,2);};
      virtual double GetE2(double x,double y){return GetData(x,y,3);};
      virtual double GetImpurity(double x,double y){return GetData(x,y,0);};

   
      virtual void SetImpurity(TF2 * Im);
      virtual void SetImpurity(double density);

      ClassDef(XY,1);

   protected:

      virtual void SOR2(int idx,bool elec); 
      double *fE2,*fC2,*fDistanceToLeft,*fDistanceToRight;//left and right are for y axis
      virtual int FindIdx(double tarx,double tary
            ,int ybegin,int yend);


      double GetData(double tarx,double tary,int thing);
      virtual void SetStepLength(double steplength1,double steplength2);
};
#endif

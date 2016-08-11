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

      virtual void SetStepLength(double steplength1,double steplength2,double steplength3);

      virtual void SaveField(const char *fout=NULL);
      virtual void LoadField(const char *fin=NULL);

      virtual double GetData(double tarx,double tary,double tarz,int thing);
            virtual void SetImpurity(TF3 * Im);
      
      virtual double GetPotential(double tarx,double tary,double tarz){return GetData(tarx,tary,tarz,1);};
      virtual double GetE1(double tarx,double tary,double tarz){return GetData(tarx,tary,tarz,2);};
      virtual double GetE2(double tarx,double tary,double tarz){return GetData(tarx,tary,tarz,3);};
      virtual double GetE3(double tarx,double tary,double tarz){return GetData(tarx,tary,tarz,4);};
      virtual double GetImpurity(double tarx,double tary,double tarz){return GetData(tarx,tary,tarz,0);};
      ClassDef(XYZ,1);

   protected:

      double *fE3,*fC3,*fDistanceToUp,*fDistanceToDown;//left and right are for y axis
      virtual int FindIdx(double tarx,double tary,
            double tarz,int begin,int end);
      virtual void SOR2(int idx,bool elec); 
};

#endif

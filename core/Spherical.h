#ifndef GEFICA_HALFBALL_H
#define GEFICA_HALFBALL_H

#include "XYZ.h"

namespace GEFICA { class Spherical; }

class GEFICA::Spherical : public GEFICA::XYZ
{
   public:

   public:
      Spherical(unsigned short n_r=0, unsigned short n_theta=0,unsigned short n_phi=0): 
         XYZ(n_r, n_theta, n_phi) {}
      virtual ~Spherical() {};

      virtual void CreateGridWithFixedStepLength(double steplength);
      virtual void RK2(int idx); 

      virtual void SaveField(const char *fout=NULL);
      virtual void LoadField(const char *fin=NULL);

      virtual double GetData(double tarx,double tary,double tarz,int thing);
      virtual void SetImpurity(TF3 * Im);

      ClassDef(Spherical,1);

   protected:
      virtual int FindIdx(double tarx,double tary,
            double tarz,int zbegin,int zend);
};

#endif

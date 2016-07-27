#ifndef GEFICA_CYLINDRICAL_H
#define GEFICA_CYLINDRICAL_H

#include "XYZ.h"
class TF3;

namespace GEFICA { class Cylindrical; }

class GEFICA::Cylindrical : public GEFICA::XYZ
{
   public:
      Cylindrical(unsigned short n1, unsigned short n2,unsigned short n3): XYZ(n1,n2,n3) {};
      virtual ~Cylindrical(){};

      virtual void CreateGridWithFixedStepLength(double steplength);
      virtual void RK2(int idx); 

      virtual void SaveField(const char *fout=NULL);
      virtual void LoadField(const char *fin=NULL);

      virtual double GetData(double tarx,double tary,double tarz,int thing);
      virtual void SetImpurity(TF3 * Im);

      ClassDef(Cylindrical,1);

   protected:
      virtual int FindIdx(double tarx,double tary,
            double tarz,int begin,int end);
};

#endif

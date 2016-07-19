#ifndef GEFICA_X_H
#define GEFICA_X_H

#include <TObject.h>
class TF1;

namespace GEFICA { class X; }

class GEFICA::X: public TObject 
{
   public:
      int x; // number of steps along the 1st axis
      int MaxIterations;
      int n;
      bool * isbegin;
      double *E1, *P,*C1,*StepNext,*StepBefore,*Impurity,Csor,E0,ER;
      double XUpSum,XDownSum,Xlimit;

   public:
      X(int nx=101);
      virtual ~X();

      virtual void Create(double steplength);
      virtual bool Iterate();

      virtual void Save(const char *fout=NULL);
      virtual void Load(const char *fin=NULL);
      virtual void SetImpurity(TF1 * Im);

      virtual double GetData(double tarx,int thing); // 1:Ex 2:f 0:Impurty

      void Initialize(detector=kPlanarX);

      ClassDef(X,1);

   protected:
      virtual void Update(int idx); 
      virtual int FindIdx(double tarx,int begin,int end);
};

#endif

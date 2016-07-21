#ifndef GEFICA_X_H
#define GEFICA_X_H

#include <TObject.h>
class TF1;

namespace GEFICA { 
   class X;
   static const double epsilon = 16*8.854187817e-12;
}

class GEFICA::X: public TObject 
{
   public:
      int n1; // number of steps along the 1st axis
      int MaxIterations;
      int n; // n = n1*n2*n3
      double Csor;
      double Precision;

      enum EMethod {
         kAnalytic,
         kRK2,
         kRK4,
      }

   public:
      X(int nx=101);
      virtual ~X();

      virtual void CreateGridWithFixedStepLength(double steplength);
      virtual bool CalculateField(EMethod method=kRK2);

      virtual void SaveField(const char *fout=NULL);
      virtual void LoadField(const char *fin=NULL);
      virtual void SetImpurity(double density);
      virtual void SetImpurity(TF1 * Im);

      virtual double GetData(double tarx,int thing); // 1:Ex 2:f 0:Impurty

      ClassDef(X,1);

   protected:
      bool * fIsFixed;
      double *fE1, *fPotential,*fC1,*fDistanceToNext,*fDistanceToPrevious,*fImpurity;

      virtual void Update(int idx); 
      virtual int FindIdx(double tarx,int begin,int end);

      virtual bool Analyic(); // Analyic calculation
      virtual bool RK2(); // 2nd-order Runge-Kutta Successive Over-Relaxation
      virtual bool RK4() {return true; } // 4th-order Runge-Kutta Successive Over-Relaxation
};

#endif

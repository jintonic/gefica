/////////////////////////////
//X                        //
//                         //
//1D field                 //
//BEGIN_HTML               
//<html>                   
//<a href="../../macro/test.C" target="_blank">example macro file</a>            
//</html>                
//END_HTML                
//../macro/test.C           
/////////////////////////////
#ifndef GEFICA_X_H
#define GEFICA_X_H

#include <TObject.h>
class TF1;

namespace GEFICA { 
   enum EMethod {
      kAnalytic,
      kRK2,
      kRK4,
   };

   class X;

   static const double epsilon = 16*8.854187817e-12;
   static const double cm =1;
   static const double volt=1;
   static const double cm3=cm*cm*cm;
}

class GEFICA::X : public TObject 
{
   public:
      int n1; // number of steps along the 1st axis
      int MaxIterations; // Iteration cuts
      int n; // n = n1*n2*n3
      double Csor; // boost Iteration speed
      double Precision;

   public:
      X(int nx=101);
      //_________
      //claim a 1D field with nx grid
      virtual ~X();

      virtual void CreateGridWithFixedStepLength(double steplength);
      virtual bool CalculateField(EMethod method=kRK2);
      //___________
      //do calculate

      virtual void SaveField(const char *fout=NULL);
      virtual void LoadField(const char *fin=NULL);
      virtual void SetImpurity(double density);
      virtual void SetImpurity(TF1 * Im);

      virtual double GetData(double tarx,int thing); 
      //__________________
      // ask thingwith number: 1:Ex 2:f 0:Impurty

      ClassDef(X,1);

   protected:
      bool * fIsFixed;
      double *fE1, *fPotential,*fC1,*fDistanceToNext,*fDistanceToPrevious,*fImpurity;

      virtual int FindIdx(double tarx,int begin,int end);

      virtual bool Analyic();
      //_________________
      // Analyic calculation
      virtual void RK2(int idx,bool elec); // 2nd-order Runge-Kutta Successive Over-Relaxation
      virtual void RK4(int idx); // 4th-order Runge-Kutta Successive Over-Relaxation
};

#endif

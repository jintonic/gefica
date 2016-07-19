#ifndef GEFICA_X_H
#define GEFICA_X_H

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TObject.h>
#include <TVectorD.h>
#include <iostream>
#include <TF1.h>

namespace GEFICA { class X; }

class GEFICA::X: public TObject 
{
   public:
      int x; // number of steps along the 1st axis
      int MaxIterations;
      int n;
      bool * isbegin;
      double *E1,*C1,*P,*StepNext,*StepBefore,*Impurity,Csor,E0,ER;
      double XUpSum,XDownSum,Xlimit;

   public:
      X(int ix) {n=ix;x=ix;};
      virtual ~X() {delete[] E1;delete [] P;delete [] C1;
         delete[] StepNext;delete[] StepBefore;delete[] isbegin;delete[] Impurity;};

      virtual void Create(double steplength);
      virtual bool Iterate();

      virtual void Save(const char *fout=NULL);
      virtual void Load(const char *fin=NULL);
      virtual void SetImpurity(TF1 * Im);

      virtual double GetData(double tarx,int thing); // 1:Ex 2:f 0:Impurty

   protected:
      virtual void Update(int idx); 
      virtual int FindIdx(double tarx,int begin,int end);

   public:
      ClassDef(X,1);
};

#endif

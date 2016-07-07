#ifndef FIELD_H
#define FIELD_H

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TObject.h>
#include <TVectorD.h>
#include <iostream>
#include <TF1.h>

namespace GEFICA {
  class Field;
}
class GEFICA::Field: public TObject 
{
  public:

    Field(int ix) {n=ix;x=ix;};
    virtual ~Field() {delete[] E1;delete [] P;delete [] C1;
      delete[] StepNext;delete[] StepBefore;delete[] isbegin;delete[] Impurity;};

    int x; // number of steps along the 1st axis
    int MaxIterations;
    int n;
    bool * isbegin;
    double *E1,*C1,*P,*StepNext,*StepBefore,*Impurity,Csor,E0,ER;
    double XUpSum,XDownSum,Xlimit;
    

    virtual void Create(double steplength);
    virtual bool Iterate();
    virtual void Update(int idx); 

    virtual void Save(const char *fout=NULL);
    virtual void Load(const char *fin=NULL);
    virtual void SetImpurity(TF1 * Im);
    
    virtual int FindIdx(double tarx,int begin,int end);
    virtual double GetData(double tarx,int thing); // 1:Ex 2:f 0:Impurty

    ClassDef(Field,1);

};


#endif


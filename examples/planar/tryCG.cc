// 1.introduction of CG:
//   https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
// 2.using CG to solve Poisson's equation:
//   http://www.cs.ucsb.edu/~gilbert/cs240a/old/cs240aSpr2011/hw2/hw2.pdf
// 3.using ROOT to minimize quadratic function defined in Ref. 2:
//   https://root.cern.ch/numerical-minimization
#include "Math/GSLMinimizer.h"
#include "Math/Functor.h"
#include <TMatrixDSparse.h> 
#include <TVectorD.h> 
#include <TSystem.h> 

const int n=20;
TMatrixD *mA;
double RosenBrock(const double *xx )
{
   //cout<<gSystem->Now().AsString()<<endl;
  //ouble aa[25]={
  //  2,-1,0,0,0,
  //  -1,2,-1,0,0,
  //  0,-1,2,-1,0,
  //  0,0,-1,2,-1,
  //  0,0,0,-1,2};
  //MatrixD mA(5,5,aa);
   double ab[n];
   for (int i=0;i<n;i++) ab[i]=0;
   ab[n-1]=1000;

   TVectorD vb(n,ab);
   //vb.Print();
   TVectorD vx1(n,xx); 
   TVectorD vx2(n,xx); 
   vx1*=(*mA);
   return 0.5*vx1*vx2-vb*vx2;
}
 
int ConjugateGradient()
{
   // Choose method upon creation between:
   // kConjugateFR, kConjugatePR, kVectorBFGS,
   // kVectorBFGS2, kSteepestDescent
   cout<<"start"<<endl;
   double aa[n*n];
   for(int i=0;i<n;i++)
   {
      for (int j=0;j<n;j++)
      {
         if(aa[i+j*n]!=-1)aa[i+j*n]=0;
         if(i==j)
         {
            aa[i+j*n]=2;
            if(i>0)aa[i+j*n-1]=-1;
            if(i<n-1)aa[i+j*n+1]=-1;
         }
      }
   }
   mA=new TMatrixD(n,n,aa);
   mA->Print();

   ROOT::Math::GSLMinimizer min(  ROOT::Math::kConjugateFR);
 
   min.SetMaxFunctionCalls(1000000);
   min.SetMaxIterations(100000);
   min.SetTolerance(0.001);
 
   ROOT::Math::Functor f(&RosenBrock,n); 
   double step[n];
   
   double variable[n] ;
   for(int i=0;i<n;i++)
   {
      step[i]=1000;
      variable[i]=500;
   }
   min.SetFunction(f);
 
   // Set the free variables to be minimized!
   for (int i=0;i<n;i++)
   {
      min.SetVariable(i,Form("%d",i),variable[i],step[i]);
   }
 //min.SetVariable(0,"x",variable[0], step[0]);
 //min.SetVariable(1,"y",variable[1], step[1]);
 //min.SetVariable(2,"z",variable[2], step[2]);
 //min.SetVariable(3,"a",variable[3], step[4]);
 //min.SetVariable(4,"b",variable[4], step[4]);
 
   min.SetPrintLevel(10);
   min.Minimize(); 
 
   const double *xs = min.X();
   cout << "Minimum: f(" << xs[0] << "," << xs[1]<<","<< xs[2] <<","<< xs[3] <<","<< xs[4] << "): " 
        << RosenBrock(xs) << endl;
 
   return 0;
}

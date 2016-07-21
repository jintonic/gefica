#include <iostream>
using namespace std;

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TVectorD.h>
#include <TF1.h>

#include "X.h"
using namespace GEFICA;

X::X(int nx) : TObject(), 
   fIsFixed(0), fE1(0), fPotential(0), fC1(0), fDistanceToNext(0), fDistanceToPrevious(0), fImpurity(0)
{ n=nx;n1=nx; }

X::~X()
{
   if (fE1) delete[] fE1;
   if (fPotential) delete[] fPotential;
   if (fC1) delete[] fC1;
   if (fDistanceToNext) delete[] fDistanceToNext;
   if (fDistanceToPrevious) delete[] fDistanceToPrevious;
   if (fIsFixed) delete[] fIsFixed;
   if (fImpurity) delete[] fImpurity;
}
bool X::Analyic()
{
  cout<<"no method can use"<<endl;
  return false; 
}

void X::CreateGridWithFixedStepLength(double steplength)
{
   Csor=1;Precision=1e-7;MaxIterations=100000;

   fE1=new double[n];
   fC1=new double[n];
   fPotential=new double[n];
   fIsFixed=new bool[n];
   fDistanceToNext=new double[n];
   fDistanceToPrevious=new double[n];
   fImpurity=new double[n];
   for (int i=n;i-->0;) {
      fIsFixed[i]=false;
      fE1[i]=0;
      fC1[i]=i*steplength;
      fPotential[i]=0;
      fDistanceToNext[i]=steplength;
      fDistanceToPrevious[i]=steplength;
      fImpurity[i]=0;
   }
}

bool X::CalculateField(EMethod method)
{
   if (method==kAnalytic) return Analyic();
   else if (method==kRK4) return RK4();
   else return RK2();
}

bool X::RK2()
{
   int cnt=0,looptime=MaxIterations;
   cout<<"aaaaa"<<endl;
   while (cnt++<looptime) {
      double XUpSum=0;
      double XDownSum=0;
      for (int i=0;i<n;i++) {
         double tmp=fPotential[i];
         Update(i);
         if(tmp>0)XDownSum+=tmp;
         else XDownSum=XDownSum-tmp;
         if(fPotential[i]-tmp>0)XUpSum=XUpSum+fPotential[i]-tmp;
         else XUpSum=XUpSum+tmp-fPotential[i];
      }
      /*if(cnt%1==0)*/cout<<cnt<<"  "<<XUpSum/XDownSum<<endl;
      if (XUpSum/XDownSum<Precision){return true;}
   }
   return false;
}

void X::Update(int idx)
{
   if (fIsFixed[idx])return;
   double density=fImpurity[idx]*1.6e-19;
   double h2=fDistanceToPrevious[idx];
   double h3=fDistanceToNext[idx];
   double tmp=density/epsilon*h2*h3+(h3*fPotential[idx-1]+h2*fPotential[idx+1])/(h2+h3);
   fPotential[idx]=Csor*(tmp-fPotential[idx])+fPotential[idx];
   fE1[idx]=(fPotential[idx+1]-fPotential[idx-1])/(h2+h3);
}

int X::FindIdx(double tarx,int begin,int end)
{
   if (begin>=end)return begin;
   int mid=(begin+end)/2;
   if(fC1[mid]>=tarx)return FindIdx(tarx,begin,mid);
   else return FindIdx(tarx,mid+1,end);
}

double X::GetData(double tarx,int thing)
{
   int idx=FindIdx(tarx,0,n-1);
   if (idx==n)
   {
      switch (thing)
      {
         case 0:return fImpurity[idx];
         case 1:return fE1[idx];
         case 2:return fPotential[idx];
      }
   }
   double ab=(tarx-fC1[idx])/fDistanceToNext[idx];
   double aa=1-ab;
   switch(thing)
   {
      case 2:return fE1[idx]*ab+fE1[idx+1]*aa;
      case 1:return fPotential[idx]*ab+fC1[idx+1]*aa;
      case 0:return fImpurity[idx]*ab+fImpurity[idx+1]*aa;
   }
   return -1;
}
void X::SaveField(const char * fout)
{
   TFile * file=new TFile(fout,"recreate","data");
   TTree * tree=new TTree("t","1D");

   TVectorD v(10);

   v[7]=(double)n1;
   v[8]=1;;
   v[9]=1;
   v[0]=(double)MaxIterations;
   v[1]=(double)n;
   v[2]=Csor;
   v.Write("v");
   bool fIsFixeds;

   double E1s,C1s,Ps,StepNexts,StepBefores,impuritys;
   tree->Branch("e1",&E1s,"e1/D"); // Electric X in x
   tree->Branch("c1",&C1s,"c1/D"); // persition in x
   tree->Branch("p",&Ps,"p/D"); // electric potential
   tree->Branch("sn",&StepNexts,"StepNext/D"); // Step length to next point in x
   tree->Branch("sb",&StepBefores,"StepBefore/D"); // Step length to before point in x
   tree->Branch("ib",&fIsFixeds,"fIsFixed/O"); // check is initial point
   tree->Branch("im",&impuritys,"impurity/D"); // Impurity
   for(int i=0;i<n;i++) {
      impuritys=fImpurity[i];
      E1s=fE1[i];
      C1s=fC1[i];
      Ps=fPotential[i];
      StepNexts=fDistanceToNext[i];
      StepBefores=fDistanceToPrevious[i];
      tree->Fill();
   }
   file->Write();
   file->Close();
   delete file;
}
void X::LoadField(const char * fin)
{
   TFile *file=new TFile(fin);
   TVectorD *v1=(TVectorD*)file->Get("v");
   double * v=v1->GetMatrixArray();
   n1		=(int)	v[7];
   MaxIterations	=(int)	v[0];
   n		=(int)	v[1];
   Csor		=	v[2];


   TChain *t =new TChain("t");
   t->Add(fin);
   bool IsFixed;
   double E1,C1,fP,fStepNext,fStepBefore,fimpurity;
   t->SetBranchAddress("e1",&E1);
   t->SetBranchAddress("c1",&C1);
   t->SetBranchAddress("p",&fP);
   t->SetBranchAddress("sn",&fStepNext);
   t->SetBranchAddress("sb",&fStepBefore);
   t->SetBranchAddress("e1",&fE1);
   t->SetBranchAddress("if",&IsFixed);
   t->SetBranchAddress("im",&fimpurity);

   fE1=new double[n];
   fC1=new double[n];
   fPotential=new double[n];
   fIsFixed=new bool[n];
   fDistanceToNext=new double[n];
   fDistanceToPrevious=new double[n];
   fImpurity=new double[n];

   for (int i=0;i<n;i++) {
      t->GetEntry(i);
      fE1[i]=E1;
      fC1[i]=C1;
      fPotential[i]=fP;
      fIsFixed[i]=fIsFixed;  
      fDistanceToNext[i]=fStepNext;
      fDistanceToPrevious[i]=fStepBefore;
      fImpurity[i]=fimpurity;
   }
   file->Close();
   delete file;
}

void X::SetImpurity(double density)
{
   for(int i=n;i-->0;) fImpurity[i]=density;
}

void X::SetImpurity(TF1 * Im)
{
   for(int i=n;i-->0;) {
      fImpurity[i]=Im->Eval((double)fC1[i]);
   }
}

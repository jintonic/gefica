#include <iostream>
using namespace std;

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TVectorD.h>
#include <TF1.h>

#include "X.h"
using namespace GEFICA;
//______________________________________________________________________________
// Create grid for 1-D field calculation.
// begin_macro
// planar1d.C           
// end_macro
ClassImp(X)
X::X(int nx) : TObject(), MaxIterations(100000), Csor(1), Precision(1e-7),
   fIsFixed(0), fE1(0), fPotential(0), fC1(0), fDistanceToNext(0), fDistanceToPrevious(0), fImpurity(0)
{ 
  //claim a 1D field with nx grisd
   n=nx;
   n1=nx; 
   floaded=false;
   if (n<10) { n=11; n1=11; }
   fE1=new double[n];
   fC1=new double[n];
   fPotential=new double[n];
   fIsFixed=new bool[n];
   fDistanceToNext=new double[n];
   fDistanceToPrevious=new double[n];
   fImpurity=new double[n];
}
      

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
  // Analyic calculation
   cout<<"no method can use"<<endl;
   return false; 
}

void X::SetStepLength(double steplength)
{
//set field step length
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


void X::SetVoltage(double anode_voltage, double cathode_voltage)
{
   fIsFixed[0]=true;
   fIsFixed[n-1]=true;
   double slope = (cathode_voltage-anode_voltage)/(n-1);
   for (int i=0; i<n; i++) {
      fPotential[i]=anode_voltage+slope*i;
   }
}

bool X::CalculateField(EMethod method)
{
  //do calculate

  floaded=true;

   if (method==kAnalytic) return Analyic();
   int cnt=0;
   while (cnt++<MaxIterations) {
      double XUpSum=0;
      double XDownSum=0;
      for (int i=1;i<n-1;i++) {
         double old=fPotential[i];
         if (method==kSOR4) SOR4(i);
         else SOR2(i,0);
         if(old>0)XDownSum+=old;
         else XDownSum-=old;
         if(fPotential[i]-old>0)XUpSum+=(fPotential[i]-old);
         else XUpSum+=-(fPotential[i]-old);
      }
      if(cnt%1000==0)
         cout<<cnt<<"  "<<XUpSum/XDownSum<<" down: "<<XDownSum<<", up: "<<XUpSum<<endl;
      if (XUpSum/XDownSum<Precision){

	  return true;
      }
   }
   return false;
}

void X::SOR2(int idx,bool elec)
{
  // 2nd-order Runge-Kutta Successive Over-Relaxation
   if (fIsFixed[idx])return ;
   double density=fImpurity[idx]*1.6e-19;
   double h2=fDistanceToPrevious[idx];
   double h3=fDistanceToNext[idx];
   double tmp=-density/epsilon*h2*h3/2+(h3*fPotential[idx-1]+h2*fPotential[idx+1])/(h2+h3);
   // over-relaxation if Csor>1
   fPotential[idx]=Csor*(tmp-fPotential[idx])+fPotential[idx];

   if(elec)fE1[idx]=(fPotential[idx+1]-fPotential[idx-1])/(h2+h3);
}

void X::SOR4(int idx)
{ 
  // 4th-order Runge-Kutta Successive Over-Relaxation
  if (fIsFixed[idx])return;

   double density=fImpurity[idx]*1.6e-19;
   double h2=fDistanceToPrevious[idx];
   double h3=fDistanceToNext[idx];
   double h1=h2;
   double xp2,xm2,xm1,xp1;
   xm1=fPotential[idx-1];
   xp1=fPotential[idx+1];
   if(idx>1)xm2=fPotential[idx-2];
   else {SOR2(idx,0);return; } 
   if(idx<n-2)xp2=fPotential[idx+2];
   else {SOR2(idx,0);return;}
   double tmp=(-1/12*xp2+4/3*xp1+4/3*xm1-1/12*xm2-density/epsilon*h1*h1)*2/5;
   fPotential[idx]=Csor*(tmp-fPotential[idx])+fPotential[idx];

   fE1[idx]=(fPotential[idx+1]-fPotential[idx-1])/(h2+h3);
}

int X::FindIdx(double tarx,int begin,int end)
{
  //search using binary search
   if (begin>=end)return begin;
   int mid=(begin+end)/2;
   if(fC1[mid]>=tarx)return FindIdx(tarx,begin,mid);
   else return FindIdx(tarx,mid+1,end);
}

double X::GetData(double tarx,int thing)
{
  // ask thingwith number: 1:Ex 2:f 0:Impurty
   int idx=FindIdx(tarx,0,n-1);
   if (idx==n)
   {
      switch (thing)
      {
         case 0:return fImpurity[idx];
         case 2:return fE1[idx];
         case 1:return fPotential[idx];
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
   tree->Branch("sb",&StepBefores,"Step)Before/D"); // Step length to before point in x
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
  floaded=true;
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
   for(int idx=0;idx-->n;)SOR2(idx,1);
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

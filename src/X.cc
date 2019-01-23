#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TVectorD.h>
#include <TStopwatch.h>
#include <TF3.h>
#include <Math/Functor.h>

#include "X.h"
#include "Units.h"
using namespace GeFiCa;
using namespace std;

X::X(int nx) : TNamed("X","X"), n1(nx), n(nx), Csor(1.95), Precision(1e-7),
   MaxIterations(100000), V0(0), V1(2000*volt)
{ 
   if (n<10) { n=11; n1=11; }

   fIsDepleted=new bool[n];
   for(int i=0;i<n;i++)fIsDepleted[i]=true;

   fIsLoaded=false;
   fE1=new double[n];
   fC1=new double[n];
   fPotential=new double[n];
   fIsFixed=new bool[n];
   fdC1p=new double[n];
   fdC1m=new double[n];
   fImpurity=new double[n];
}
//_____________________________________________________________________________
//
X::~X()
{
   if (fE1) delete[] fE1;
   if (fPotential) delete[] fPotential;
   if (fC1) delete[] fC1;
   if (fdC1p) delete[] fdC1p;
   if (fdC1m) delete[] fdC1m;
   if (fIsFixed) delete[] fIsFixed;
   if (fImpurity) delete[] fImpurity;
   if (fIsDepleted) delete[] fIsDepleted;
}
//_____________________________________________________________________________
//
bool X::Analytic()
{
   Printf("There is no analytic solution for this setup");
   return false; 
}
//_____________________________________________________________________________
//
//X& X::operator=(GeFiCa::X *anotherfield)
//{
//    return (GeFiCa::X)anotherfield.Clone("newfield");
//}
//_____________________________________________________________________________
//
X& X::operator+=(GeFiCa::X *anotherfield)
{
   if (n!=anotherfield->n) {
      Warning("Add", 
            "Only same type of detector can be added together! Do nothing.");
      return *this; 
   }
   for (int i=0; i<n; i++)
   {
      fPotential[i]=fPotential[i]+anotherfield->fPotential[i];
      fImpurity[i]+=anotherfield->fImpurity[i];
   }
   V0+=anotherfield->V0; V1+=anotherfield->V1; 
   return *this;
}
//_____________________________________________________________________________
//
X& X::operator*=(double p)
{
   for (int i=0; i<n; i++) fPotential[i]=fPotential[i]*p;
   V0*=p; V1*=p;
   return *this;
}
//_____________________________________________________________________________
//
int X::Findmax()
{
   double max=fPotential[0];
   int maxn=0;
   for(int i=1;i<n;i++)
   {
      if(fPotential[i]>max)
      {
         maxn=i;
         max=fPotential[i];
      }
   }
   return maxn;
}
//_____________________________________________________________________________
//
int X::Findmin()
{
   double min=fPotential[0];
   int minn=0;
   for(int i=1;i<n;i++)
   {
      if(fPotential[i]<min)
      {
         minn=i;
         min=fPotential[i];
      }
   }
   return minn;
}
//_____________________________________________________________________________
//
bool X::IsDepleted()
{
  for(int i=0;i<n;i++)
  {
     SOR2(i,0);
     if (!fIsDepleted[i])return false;
  }
  return true;
}
//_____________________________________________________________________________
//
void X::SetStepLength(double steplength)
{
   //set field step length
   for (int i=n;i-->0;) {
      fIsFixed[i]=false;
      fC1[i]=i*steplength;
      fdC1p[i]=steplength;
      fdC1m[i]=steplength;
   }
}
//_____________________________________________________________________________
//
int* X::FindSurroundingMatrix(int idx)
{
   int *tmp=new int[3];
   tmp[0]=idx;
   if(idx-1<0)tmp[1]=1;
   else tmp[1]=idx-1;
   if(idx+1>=n)tmp[2]=n-2;
   else tmp[2]=idx+1;
   return tmp;
}
//_____________________________________________________________________________
//
bool X::CalculatePotential(EMethod method)
{
   TStopwatch watch; watch.Start();
   if (method==kAnalytic) return Analytic();
   cout<<" Calculate field ..."<<endl;
   int cnt=0;
   bool NotImpurityPotential;
   if (V0==0&&V1==0)NotImpurityPotential=false;
   else NotImpurityPotential=true;
   while (cnt++<MaxIterations) {
      double XUpSum=0;
      double XDownSum=0;
      for (int i=n-1;i>=0;i--) {
         double old=fPotential[i];
         SOR2(i,NotImpurityPotential);
         if(old>0)XDownSum+=old;
         else XDownSum-=old;
         double diff=fPotential[i]-old;
         if(diff>0)XUpSum+=(diff);
         else XUpSum-=(diff);
      }
      double cp = XUpSum/XDownSum; // current precision
      if (cnt%100==0)
         Printf("  %05d, target precision: %e, current precision: %e", 
               cnt, Precision, cp);

      if (cp<Precision) {
         for (int i=0; i<n; i++) if (!CalculateField(i)) return false;
         Printf("  %05d, target precision: %e, current precision: %e", 
               cnt, Precision, cp);
         cout<<" Done. Spent "; watch.Stop(); watch.Print();
         return true;
      }
   }
   return false;
}
//_____________________________________________________________________________
//
void X::SOR2(int idx,bool NotImpurityPotential)
{
   // 2nd-order Runge-Kutta Successive Over-Relaxation
   if (fIsFixed[idx])return ;
   double density=-fImpurity[idx]*Qe;
   double h2=fdC1m[idx];
   double h3=fdC1p[idx];
   double p2=fPotential[idx-1];
   double p3=fPotential[idx+1];

   double tmp=-density/epsilon*h2*h3/2+(h3*fPotential[idx-1]+h2*fPotential[idx+1])/(h2+h3);
   
   //find minmium and maxnium of all five grid, the new one should not go overthem.
   //find min
   double min=p2;
   double max=p2;
   if(min>p3)min=p3;
   //find max
   if(max<p3)max=p3;
//if tmp is greater or smaller than max and min, set tmp to it.
   
   //fPotential[idx]=Csor*(tmp-fPotential[idx])+fPotential[idx];
   double oldP=fPotential[idx];
      tmp=Csor*(tmp-oldP)+oldP;
   
   if(tmp<min)
   {
      fPotential[idx]=min;
      fIsDepleted[idx]=false;
   }
   else if(tmp>max)
   {
      fPotential[idx]=max;
      fIsDepleted[idx]=false;
   }
   else
      fIsDepleted[idx]=true;
   if(fIsDepleted[idx]||!NotImpurityPotential)
   {
      //over relax
      fPotential[idx]=tmp;
   }
}
//_____________________________________________________________________________
//
int X::FindIdx(double tarx,int begin,int end)
{
   //search using binary search
   if (begin>=end)return end;
   int mid=(begin+end)/2;
   if(fC1[mid]>=tarx)return FindIdx(tarx,begin,mid);
   else return FindIdx(tarx,mid+1,end);
}

//_____________________________________________________________________________
//
double X::GetData(double tarx, EOutput output)
{
   int idx=FindIdx(tarx,0,n-1);
   if (idx==n) {
      switch (output) {
         case 0:return fImpurity[idx];
         case 2:return fE1[idx];
         case 1:return fPotential[idx];
         default: return -1;
      }
   }
   double ab=(-tarx+fC1[idx])/fdC1p[idx];
   double aa=1-ab;
   switch(output) {
      case 2:return fE1[idx]*ab+fE1[idx-1]*aa;
      case 1:return fPotential[idx]*ab+fC1[idx-1]*aa;
      case 0:return fImpurity[idx]*ab+fImpurity[idx-1]*aa;
      default: return -1;
   }
   return -1;
}
//_____________________________________________________________________________
//
void X::SaveField(const char * fout)
{
   TFile * file=new TFile(fout,"recreate","data");
   TTree * tree=new TTree("t","1D");

   TVectorD v(10);

   v[7]=(double)n1;
   v[8]=1;
   v[9]=1;
   v[0]=(double)MaxIterations;
   v[1]=(double)n;
   v[2]=Csor;
   v.Write("v");
   bool fIsFixeds;

   double E1s,C1s,Ps,dC1p,dC1m,impuritys;
   tree->Branch("e1",&E1s,"e1/D"); // Electric X in x
   tree->Branch("c1",&C1s,"c1/D"); // persition in x
   tree->Branch("p",&Ps,"p/D"); // electric potential
   tree->Branch("dC1p",&dC1p,"dC1p/D"); // Step length to next point in x
   tree->Branch("dC1m",&dC1m,"dC1m/D"); // Step length to before point in x
   tree->Branch("ib",&fIsFixeds,"ib/O"); // check is initial point
   tree->Branch("im",&impuritys,"im/D"); // Impurity
   for(int i=0;i<n;i++) {
      impuritys=fImpurity[i];
      E1s=fE1[i];
      C1s=fC1[i];
      Ps=fPotential[i];
      dC1p=fdC1p[i];
      dC1m=fdC1m[i];
      fIsFixeds=fIsFixed[i];
      tree->Fill();
   }
   file->Write();
   file->Close();
   delete file;
}
//_____________________________________________________________________________
//
void X::LoadField(const char * fin)
{
   //will calculate electric field after load
   fIsLoaded=true;
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
   double E1,C1,P,dC1p,dC1m,fimpurity;
   t->SetBranchAddress("e1",&E1);
   t->SetBranchAddress("c1",&C1);
   t->SetBranchAddress("p",&P);
   t->SetBranchAddress("dC1p",&dC1p);
   t->SetBranchAddress("dC1m",&dC1m);
   t->SetBranchAddress("ib",&IsFixed);
   t->SetBranchAddress("im",&fimpurity);

   fE1=new double[n];
   fC1=new double[n];
   fPotential=new double[n];
   fIsFixed=new bool[n];
   fdC1p=new double[n];
   fdC1m=new double[n];
   fImpurity=new double[n];

   for (int i=0;i<n;i++) {
      t->GetEntry(i);
      fE1[i]=E1;
      fC1[i]=C1;
      fPotential[i]=P;
      fIsFixed[i]=IsFixed;  
      fdC1p[i]=dC1p;
      fdC1m[i]=dC1m;
      fImpurity[i]=fimpurity;
   }
   file->Close();
   delete file;
   for(int idx=0;idx-->n;)SOR2(idx,1);
}
//_____________________________________________________________________________
//
void X::SetImpurity(TF3 *fi1)
{
   Initialize();
   for (int i=n;i-->0;) fImpurity[i]=fi1->Eval(fC1[i]);
}
//_____________________________________________________________________________

bool X::CalculateField(int idx)
{
   if (fdC1p[idx]==0 || fdC1m[idx]==0) return false;

   if (idx%n1==0) // C1 lower boundary
      fE1[idx]=(fPotential[idx]-fPotential[idx+1])/fdC1p[idx];
   else if (idx%n1==n1-1) // C1 upper boundary
      fE1[idx]=(fPotential[idx]-fPotential[idx-1])/fdC1m[idx];
   else { // bulk
      fE1[idx]=(fPotential[idx-1]-fPotential[idx+1])/(fdC1m[idx]+fdC1p[idx]);
   }
   return true;
}

//_____________________________________________________________________________
//
double X::GetCapacitance()
{
   cout<<"Calculate detector capacitance..."<<endl;
   // set impurity to zero
   double *tmpImpurity=fImpurity;
   for (int i=0;i<n;i++) {
      if (fImpurity[i]!=0) {
         //impurity not clear,return
         //return -1;
         fImpurity=new double[n];
         for (int j=0;j<n;j++) {
            fImpurity[j]=0;
            if (!fIsFixed[j] && !fIsDepleted[j]) fIsFixed[j]=true;
         }
         break;
      }
   }
   // calculate potential without impurity
   CalculatePotential(GeFiCa::kSOR2);
   // set impurity back
   if(fImpurity!=tmpImpurity) delete []fImpurity;
   fImpurity=tmpImpurity;

   // calculate C based on CV^2/2 = epsilon int E^2 dx^3 / 2
   double dV=V0-V1; if(dV<0)dV=-dV;
   double SumofElectricField=0;
   for(int i=0;i<n;i++) {
      SumofElectricField+=fE1[i]*fE1[i]*fdC1p[i]*cm*cm;
      if (!fIsDepleted[i]) fIsFixed[i]=false;
   }
   cout<<"Done."<<endl;
   return SumofElectricField*epsilon/dV/dV;
}

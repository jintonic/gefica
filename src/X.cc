#include <TF3.h>
#include <TTree.h>
#include <TStyle.h>
#include <TStopwatch.h>

#include "X.h"
#include "Units.h"
using namespace GeFiCa;

X::X(int nx, const char *name, const char *title) : TNamed(name,title), V0(0),
   V1(2e3*volt), n1(nx), n(nx), MaxIterations(5000), Nsor(0), Csor(1.95),
   Precision(1e-7*volt)
{
   if (n<10) { Warning("X","n<10, set it to 11"); n=11; n1=11; }

   fV=new double[n];
   fE1=new double[n]; fE2=new double[n]; fE3=new double[n];
   fC1=new double[n]; fC2=new double[n]; fC3=new double[n];

   fdC1p=new double[n]; fdC1m=new double[n];
   fdC2p=new double[n]; fdC2m=new double[n];
   fdC3p=new double[n]; fdC3m=new double[n];

   fIsFixed=new bool[n];
   fIsDepleted=new bool[n];
   fImpurity=new double[n];

   for (int i=0;i<n;i++) {
      fV[i]=0;
      fE1[i]=0; fE2[i]=0; fE3[i]=0;
      fC1[i]=0; fC2[i]=0; fC3[i]=0;

      fdC1m[i]=0; fdC1p[i]=0;
      fdC2m[i]=0; fdC2p[i]=0;
      fdC3m[i]=0; fdC3p[i]=0;

      fIsFixed[i]=false;
      fIsDepleted[i]=true;
      fImpurity[i]=0;
   }

   fTree=NULL; fImpDist=NULL;

   gROOT->SetStyle("Plain"); // pick up a good default style to modify
   gStyle->SetLegendBorderSize(0);
   gStyle->SetLegendFont(132);
   gStyle->SetLabelFont(132,"XYZ");
   gStyle->SetTitleFont(132,"XYZ");
   gStyle->SetLabelSize(0.05,"XYZ");
   gStyle->SetTitleSize(0.05,"XYZ");
   gStyle->SetTitleOffset(-0.4,"Z");
   gStyle->SetPadTopMargin(0.02);
}
//_____________________________________________________________________________
//
X::~X()
{
   if (fV) delete[] fV;
   if (fE1) delete[] fE1;
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
   Info("Analytic", "There is no analytic solution for this setup");
   return false; 
}
//_____________________________________________________________________________
//
X& X::operator+=(GeFiCa::X *other)
{
   if (n!=other->n) {
      Warning("+=", 
            "Only same type of detector can be added together! Do nothing.");
      return *this; 
   }
   for (int i=0; i<n; i++) {
      fV[i]=fV[i]+other->fV[i];
      fImpurity[i]+=other->fImpurity[i];
   }
   V0+=other->V0; V1+=other->V1; 
   return *this;
}
//_____________________________________________________________________________
//
X& X::operator*=(double p)
{
   for (int i=0; i<n; i++) fV[i]=fV[i]*p;
   V0*=p; V1*=p;
   return *this;
}
//_____________________________________________________________________________
//
int X::GetIdxOfMaxV()
{
   double max=fV[0];
   int maxn=0;
   for(int i=1;i<n;i++) {
      if(fV[i]>max) {
         maxn=i;
         max=fV[i];
      }
   }
   return maxn;
}
//_____________________________________________________________________________
//
int X::GetIdxOfMinV()
{
   double min=fV[0];
   int minn=0;
   for(int i=1;i<n;i++) {
      if(fV[i]<min) {
         minn=i;
         min=fV[i];
      }
   }
   return minn;
}
//_____________________________________________________________________________
//
bool X::IsDepleted()
{
   for(int i=0;i<n;i++) {
      DoSOR2(i); // calculate one more time in case of 
      //adding two fields together, one is depleted, the other is not
      if (!fIsDepleted[i]) return false;
   }
   return true;
}
//_____________________________________________________________________________
//
void X::SetStepLength(double stepLength)
{
   for (int i=n;i-->0;) {
      fIsFixed[i]=false;
      fC1[i]=i*stepLength;
      fdC1p[i]=stepLength;
      fdC1m[i]=stepLength;
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
   if (fdC1p[0]==0) Initialize(); // setup and initialize grid if it's not done
   if (fImpDist && fImpurity[0]==0) // set impurity values if it's not done yet
      for (int i=n;i-->0;) fImpurity[i]=fImpDist->Eval(fC1[i], fC2[i], fC3[i]);
   if (method==kAnalytic) return Analytic();

   Info("CalculatePotential","Start SOR...");
   TStopwatch watch; watch.Start();
   double cp=1; // current presision
   while (Nsor<MaxIterations) {
      if (Nsor%100==0) Printf("%4d steps, precision: %.1e (target: %.0e)", 
               Nsor, cp, Precision);
      double XUpSum=0;
      double XDownSum=0;
      for (int i=n-1;i>=0;i--) {
         double old=fV[i];
         DoSOR2(i);
         if(old>0)XDownSum+=old;
         else XDownSum-=old;
         double diff=fV[i]-old;
         if(diff>0)XUpSum+=(diff);
         else XUpSum-=(diff);
      }
      cp = XUpSum/XDownSum;
      Nsor++;
      if (cp<Precision) break;
   }
   for (int i=0; i<n; i++) if (!CalculateField(i)) return false;
   Printf("%4d steps, precision: %.1e (target: %.0e)", Nsor, cp, Precision);
   Info("CalculatePotential", "CPU time: %.1f s", watch.CpuTime());
   return true;
}
//_____________________________________________________________________________
//
void X::DoSOR2(int idx)
{
   // 2nd-order Runge-Kutta Successive Over-Relaxation
   if (fIsFixed[idx])return ;
   double rho=-fImpurity[idx]*Qe;
   double h2=fdC1m[idx];
   double h3=fdC1p[idx];
   double p2=fV[idx-1];
   double p3=fV[idx+1];

   double tmp=-rho/epsilon*h2*h3/2
      + (h3*fV[idx-1]+h2*fV[idx+1])/(h2+h3);

   //find minmium and maxnium of all five grid, the new one should not go overthem.
   //find min
   double min=p2;
   double max=p2;
   if(min>p3)min=p3;
   //find max
   if(max<p3)max=p3;
   //if tmp is greater or smaller than max and min, set tmp to it.

   //fV[idx]=Csor*(tmp-fV[idx])+fV[idx];
   double oldP=fV[idx];
   tmp=Csor*(tmp-oldP)+oldP;

   if(tmp<min) {
      fV[idx]=min;
      fIsDepleted[idx]=false;
   } else if(tmp>max) {
      fV[idx]=max;
      fIsDepleted[idx]=false;
   } else
      fIsDepleted[idx]=true;

   if(fIsDepleted[idx]||V0==V1) fV[idx]=tmp;
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
         case 1:return fV[idx];
         default: return -1;
      }
   }
   double ab=(-tarx+fC1[idx])/fdC1p[idx];
   double aa=1-ab;
   switch(output) {
      case 2:return fE1[idx]*ab+fE1[idx-1]*aa;
      case 1:return fV[idx]*ab+fC1[idx-1]*aa;
      case 0:return fImpurity[idx]*ab+fImpurity[idx-1]*aa;
      default: return -1;
   }
   return -1;
}
//_____________________________________________________________________________
//
bool X::CalculateField(int idx)
{
   if (fdC1p[idx]==0 || fdC1m[idx]==0) return false;

   if (idx%n1==0) // C1 lower boundary
      fE1[idx]=(fV[idx]-fV[idx+1])/fdC1p[idx];
   else if (idx%n1==n1-1) // C1 upper boundary
      fE1[idx]=(fV[idx]-fV[idx-1])/fdC1m[idx];
   else // bulk
      fE1[idx]=(fV[idx-1]-fV[idx+1])/(fdC1m[idx]+fdC1p[idx]);

   return true;
}
//_____________________________________________________________________________
//
double X::GetCapacitance()
{
   Info("GetCapacitance","Start...");
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
   Info("GetCapacitance","Done.");
   return SumofElectricField*epsilon/dV/dV;
}
//_____________________________________________________________________________
//
TTree* X::GetTree(bool createNew)
{
   if (fTree) { if (createNew) delete fTree; else return fTree; }

   // define tree
   bool b,d; double v,te,e1,e2,e3,c1,c2,c3;
   fTree = new TTree("t","field data");
   fTree->Branch("potential",&v,"v/D");
   fTree->Branch("total E  ",&te,"e/D");
   // 1D data
   fTree->Branch("E1            ",&e1,"e1/D");
   fTree->Branch("1st coordinate",&c1,"c1/D");
   // initialize values
   if (fdC1p[0]==0) Initialize(); // setup & initialize grid
   if (fImpDist && fImpurity[0]==0) for (int i=n;i-->0;)
      fImpurity[i]=fImpDist->Eval(fC1[i]/cm, fC2[i]/cm, fC3[i]/cm);

   if (fdC2p[0]!=0) { // if it is a 2D grid
      fTree->Branch("E2            ",&e2,"e2/D");
      fTree->Branch("2nd coordinate",&c2,"c2/D");
   }
   if (fdC3p[0]!=0) { // if it is a 3D grid
      fTree->Branch("E3            ",&e3,"e3/D");
      fTree->Branch("3rd coordinate",&c3,"c3/D");
   }
   fTree->Branch("boundary flag",&b,"b/O"); // boundary flag
   fTree->Branch("depletion flag",&d,"d/O"); // depletion flag

   // fill tree
   Info("GetTree","%d entries",n);
   for (int i=0; i<n; i++) {
      e1= fE1[i]; c1= fC1[i]; // 1D data
      if (fdC2p[i]!=0) { e2=fE2[i]; c2=fC2[i]; } // 2D data
      if (fdC3p[i]!=0) { e3=fE3[i]; c3=fC3[i]; } // 3D data
      v = fV[i]; b = fIsFixed[i]; d = fIsDepleted[i]; // common data
      if (fdC3p[i]!=0) te=TMath::Sqrt(e1*e1 + e2*e2 + e3*e3);
      else { if (fdC2p[i]!=0) te=TMath::Sqrt(e1*e1+e2*e2); else te=e1; }
      fTree->Fill();
   }

   fTree->GetListOfBranches()->ls();
   gDirectory->ls();
   fTree->ResetBranchAddresses(); // disconnect from local variables
   return fTree;
}

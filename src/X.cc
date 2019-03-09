#include <TF3.h>
#include <TTree.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TStopwatch.h>

#include "X.h"
#include "Units.h"
using namespace GeFiCa;

X::X(int nx, const char *name, const char *title) : TNamed(name,title), V0(0),
   V1(2e3*volt), MaxIterations(5000), Csor(1.94), Precision(1e-7*volt),
   Gsor(0), fN(nx), fN1(nx), fN2(0), fN3(0), fTree(0), fImpDist(0)
{
   if (fN<10) { Warning("X","fN<10, set it to 11"); fN=11; fN1=11; }

   fV=new double[fN];
   fE1=new double[fN]; fE2=new double[fN]; fE3=new double[fN];
   fC1=new double[fN]; fC2=new double[fN]; fC3=new double[fN];

   fdC1p=new double[fN]; fdC1m=new double[fN];
   fdC2p=new double[fN]; fdC2m=new double[fN];
   fdC3p=new double[fN]; fdC3m=new double[fN];

   fIsFixed=new bool[fN];
   fIsDepleted=new bool[fN];
   fImpurity=new double[fN];

   for (int i=0;i<fN;i++) {
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

   // pick up a good style to modify
   gROOT->SetStyle("Plain");
   gStyle->SetLegendBorderSize(0);
   gStyle->SetLegendFont(132);
   gStyle->SetLabelFont(132,"XYZ");
   gStyle->SetTitleFont(132,"XYZ");
   gStyle->SetLabelSize(0.05,"XYZ");
   gStyle->SetTitleSize(0.05,"XYZ");
   gStyle->SetTitleOffset(-0.4,"Z");
   gStyle->SetPadTopMargin(0.02);
   // create a smoother palette than the default one
   const int nRGBs = 5;
   const int nCont = 255;
   double stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
   double red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
   double green[nRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
   double blue[nRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
   TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
   gStyle->SetNumberContours(nCont);
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
   if (fN!=other->fN) {
      Warning("+=", 
            "Only same type of detector can be added together! Do nothing.");
      return *this; 
   }
   for (int i=0; i<fN; i++) {
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
   for (int i=0; i<fN; i++) fV[i]=fV[i]*p;
   V0*=p; V1*=p;
   return *this;
}
//_____________________________________________________________________________
//
int X::GetIdxOfMaxV()
{
   double max=fV[0];
   int maxn=0;
   for(int i=1;i<fN;i++) {
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
   for(int i=1;i<fN;i++) {
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
   for(int i=0;i<fN;i++) {
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
   for (int i=fN;i-->0;) {
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
   if(idx+1>=fN)tmp[2]=fN-2;
   else tmp[2]=idx+1;
   return tmp;
}
//_____________________________________________________________________________
//
bool X::CalculatePotential(EMethod method)
{
   if (fdC1p[0]==0) Initialize(); // setup and initialize grid if it's not done
   if (method==kAnalytic) return Analytic();

   Info("CalculatePotential","Start SOR...");
   if (Gsor==0) {
      Gsor = new TGraph; Gsor->SetName("Gsor");
      Gsor->SetTitle(";Number of iterations;log10(precision)");
   }
   else Gsor->Set(0); // reset the graph
   double cp=1; // current presision
   int it=0; // # of iterations
   TStopwatch watch; watch.Start();
   while (it<MaxIterations) {
      if (it%100==0) {
         Printf("%4d steps, precision: %.1e (target: %.0e)", 
               it, cp, Precision);
         if (it!=0) Gsor->SetPoint(Gsor->GetN(),it,TMath::Log10(cp));
      }
      double XUpSum=0;
      double XDownSum=0;
      for (int i=fN-1;i>=0;i--) {
         double old=fV[i];
         DoSOR2(i);
         if(old>0)XDownSum+=old;
         else XDownSum-=old;
         double diff=fV[i]-old;
         if(diff>0)XUpSum+=(diff);
         else XUpSum-=(diff);
      }
      cp = XUpSum/XDownSum;
      it++;
      if (cp<Precision) break;
   }
   for (int i=0; i<fN; i++) if (!CalculateField(i)) return false;
   Printf("%4d steps, precision: %.1e (target: %.0e)", it, cp, Precision);
   Gsor->SetPoint(Gsor->GetN(),it,TMath::Log10(cp));
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
double X::GetData(double x, double y, double z, double *data)
{
   int idx=FindIdx(x,0,fN-1);
   if (idx==fN) return data[idx];
   double ab=(-x+fC1[idx])/fdC1p[idx];
   double aa=1-ab;
   return data[idx]*ab+data[idx-1]*aa;
}
//_____________________________________________________________________________
//
bool X::CalculateField(int idx)
{
   if (fdC1p[idx]==0 || fdC1m[idx]==0) return false;

   if (idx%fN1==0) // C1 lower boundary
      fE1[idx]=(fV[idx]-fV[idx+1])/fdC1p[idx];
   else if (idx%fN1==fN1-1) // C1 upper boundary
      fE1[idx]=(fV[idx]-fV[idx-1])/fdC1m[idx];
   else // bulk
      fE1[idx]=(fV[idx-1]-fV[idx+1])/(fdC1m[idx]+fdC1p[idx]);

   return true;
}
//_____________________________________________________________________________
//
double X::GetC()
{
   Info("GetC","Start...");
   CalculatePotential(GeFiCa::kSOR2); // identify undepleted region
   // set impurity to zero
   double *tmpImpurity=fImpurity;
   for (int i=0;i<fN;i++) {
      if (fImpurity[i]!=0) {
         fImpurity=new double[fN];
         for (int j=0;j<fN;j++) {
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
   for(int i=0;i<fN;i++) {
      SumofElectricField+=fE1[i]*fE1[i]*fdC1p[i]*cm*cm;
      if (!fIsDepleted[i]) fIsFixed[i]=false;
   }
   double c=SumofElectricField*epsilon/dV/dV;
   Info("GetC","%.2f pF",c/pF);
   return c;
}
//_____________________________________________________________________________
//
TTree* X::GetTree(bool createNew)
{
   if (fTree) { if (createNew) delete fTree; else return fTree; }

   // define tree
   bool b,d; double v,im,te,e1,e2,e3,c1,c2,c3;
   fTree = new TTree("t","field data");
   fTree->Branch("potential",&v,"v/D");
   fTree->Branch("total E  ",&te,"e/D");
   // 1D data
   fTree->Branch("E1            ",&e1,"e1/D");
   fTree->Branch("1st coordinate",&c1,"c1/D");
   // initialize values
   if (fdC1p[0]==0) Initialize(); // setup & initialize grid

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
   fTree->Branch("impurity",&im,"im/D"); // depletion flag

   // fill tree
   Info("GetTree","%d entries",fN);
   for (int i=0; i<fN; i++) {
      e1= fE1[i]; c1= fC1[i]; // 1D data
      if (fdC2p[i]!=0) { e2=fE2[i]; c2=fC2[i]; } // 2D data
      if (fdC3p[i]!=0) { e3=fE3[i]; c3=fC3[i]; } // 3D data
      v = fV[i]; b = fIsFixed[i]; d = fIsDepleted[i]; im=fImpurity[i]; // common data
      if (fdC3p[i]!=0) te=TMath::Sqrt(e1*e1 + e2*e2 + e3*e3);
      else { if (fdC2p[i]!=0) te=TMath::Sqrt(e1*e1+e2*e2); else te=e1; }
      fTree->Fill();
   }

   fTree->GetListOfBranches()->ls();
   gDirectory->ls();
   fTree->ResetBranchAddresses(); // disconnect from local variables
   return fTree;
}
//_____________________________________________________________________________
//
void X::SetGridImpurity()
{
   if (fImpDist && fImpurity[0]==0) // set impurity values if it's not done yet
      for (int i=fN;i-->0;) fImpurity[i]=fImpDist->Eval(fC1[i], fC2[i], fC3[i]);
}
//_____________________________________________________________________________
//
int X::GetNsor()
{
   if (Gsor) return Gsor->GetX()[Gsor->GetN()-1];
   else return 0;
}

#include <TTree.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TStopwatch.h>

#include "Grid.h"
#include "Units.h"
using namespace GeFiCa;

X::X(int nx, const char *name, const char *title) : TNamed(name,title), Bias[0](0),
   Bias[1](2e3*volt), MaxIterations(5000), RelaxationFactor(1.94), Precision(1e-7*volt),
   Gsor(0), fN(nx), fN1(nx), fN2(0), fN3(0), fTree(0), fImpDist(0)
{
   if (fN<10) { Warning("X","fN<10, set it to 11"); fN=11; fN1=11; }

   V=new double[fN];
   E1=new double[fN]; E2=new double[fN]; E3=new double[fN];
   C1=new double[fN]; C2=new double[fN]; C3=new double[fN];

   dC1p=new double[fN]; dC1m=new double[fN];
   dC2p=new double[fN]; dC2m=new double[fN];
   dC3p=new double[fN]; dC3m=new double[fN];

   fIsFixed=new bool[fN];
   fIsDepleted=new bool[fN];
   fImpurity=new double[fN];

   for (int i=0;i<fN;i++) {
      V[i]=0;
      E1[i]=0; E2[i]=0; E3[i]=0;
      C1[i]=0; C2[i]=0; C3[i]=0;

      dC1m[i]=0; dC1p[i]=0;
      dC2m[i]=0; dC2p[i]=0;
      dC3m[i]=0; dC3p[i]=0;

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
   if (V) delete[] V;
   if (E1) delete[] E1;
   if (C1) delete[] C1;
   if (dC1p) delete[] dC1p;
   if (dC1m) delete[] dC1m;
   if (fIsFixed) delete[] fIsFixed;
   if (fImpurity) delete[] fImpurity;
   if (fIsDepleted) delete[] fIsDepleted;
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
      V[i]=V[i]+other->V[i];
      fImpurity[i]+=other->fImpurity[i];
   }
   Bias[0]+=other->Bias[0]; Bias[1]+=other->Bias[1]; 
   return *this;
}
//_____________________________________________________________________________
//
X& X::operator*=(double p)
{
   for (int i=0; i<fN; i++) V[i]=V[i]*p;
   Bias[0]*=p; Bias[1]*=p;
   return *this;
}
//_____________________________________________________________________________
//
int X::GetIdxOfMaxV()
{
   double max=V[0];
   int maxn=0;
   for(int i=1;i<fN;i++) {
      if(V[i]>max) {
         maxn=i;
         max=V[i];
      }
   }
   return maxn;
}
//_____________________________________________________________________________
//
int X::GetIdxOfMinV()
{
   double min=V[0];
   int minn=0;
   for(int i=1;i<fN;i++) {
      if(V[i]<min) {
         minn=i;
         min=V[i];
      }
   }
   return minn;
}
//_____________________________________________________________________________
//
bool X::IsDepleted()
{
   for(int i=0;i<fN;i++) {
      OverRelaxAt(i); // calculate one more time in case of 
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
      C1[i]=i*stepLength;
      dC1p[i]=stepLength;
      dC1m[i]=stepLength;
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
bool X::SuccessiveOverRelax()
{
   if (dC1p[0]==0) Initialize(); // setup and initialize grid if it's not done

   Info("SuccessiveOverRelax","Start...");
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
      for (int i=0;i<fN;i++) {
         double old=V[i];
         OverRelaxAt(i);
         if(old>0)XDownSum+=old;
         else XDownSum-=old;
         double diff=V[i]-old;
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
   Info("SuccessiveOverRelax", "CPU time: %.1f s", watch.CpuTime());
   return true;
}
//_____________________________________________________________________________
//
void X::OverRelaxAt(int idx)
{
   // 2nd-order Runge-Kutta Successive Over-Relaxation
   if (fIsFixed[idx])return ;
   double rho=-fImpurity[idx]*Qe;
   double h2=dC1m[idx];
   double h3=dC1p[idx];
   double p2=V[idx-1];
   double p3=V[idx+1];

   double tmp=-rho/epsilon*h2*h3/2
      + (h3*V[idx-1]+h2*V[idx+1])/(h2+h3);

   //find minmium and maxnium of all five grid, the new one should not go overthem.
   //find min
   double min=p2;
   double max=p2;
   if(min>p3)min=p3;
   //find max
   if(max<p3)max=p3;
   //if tmp is greater or smaller than max and min, set tmp to it.

   //V[idx]=RelaxationFactor*(tmp-V[idx])+V[idx];
   double oldP=V[idx];
   tmp=RelaxationFactor*(tmp-oldP)+oldP;

   if(tmp<min) {
      V[idx]=min;
      fIsDepleted[idx]=false;
   } else if(tmp>max) {
      V[idx]=max;
      fIsDepleted[idx]=false;
   } else
      fIsDepleted[idx]=true;

   if(fIsDepleted[idx]||Bias[0]==Bias[1]) V[idx]=tmp;
}
//_____________________________________________________________________________
//
int X::FindIdx(double tarx,int begin,int end)
{
   if (end==-1) end=fN1-1;
   //search using binary search
   if (begin>=end)return end;
   int mid=(begin+end)/2;
   if(C1[mid]>=tarx)return FindIdx(tarx,begin,mid);
   else return FindIdx(tarx,mid+1,end);
}
//_____________________________________________________________________________
//
int X::FindIdx(double tarx,double tary ,int begin,int end)
{
   if (end==-1) end=fN2-1;
   //search using binary search
   // if(begin>=end)cout<<"to x"<<begin<<" "<<end<<endl;;
   if(begin>=end)return FindIdx(tarx,end*fN1,(end+1)*fN1-1);
   int mid=((begin+end)/2);
   if(C2[mid*fN1]>=tary){//cout<<"firsthalf"<<begin<<" "<<end<<endl; 
      return FindIdx(tarx,tary,begin,mid);
   }
   else{//cout<<"senondhalf"<<begin<<" "<<end<<endl; 
      return FindIdx(tarx,tary,mid+1,end);}
}
//_____________________________________________________________________________
//
int X::FindIdx(double tarx, double tary,double tarz,int begin,int end)
{
   if (end==-1) end=fN3-1;
   //search using binary search
   if(begin>=end)return FindIdx(tarx,tary,begin,begin+fN1*fN2-1);
   int mid=((begin/(fN1*fN2)+end/(fN1*fN2))/2)*fN1*fN2;
   if(C3[mid]>=tarz)return FindIdx(tarx,tary,tarz,begin,mid);
   else return FindIdx(tarx,tary,tarz,mid+1,end);
}
//_____________________________________________________________________________
//
double X::GetData(double x, double y, double z, double *data)
{
   int idx=FindIdx(x);
   if (idx==fN) return data[idx];
   double ab=(-x+C1[idx])/dC1p[idx];
   double aa=1-ab;
   return data[idx]*ab+data[idx-1]*aa;
}
//_____________________________________________________________________________
//
bool X::CalculateField(int idx)
{
   if (dC1p[idx]==0 || dC1m[idx]==0) return false;

   if (idx%fN1==0) // C1 lower boundary
      E1[idx]=(V[idx]-V[idx+1])/dC1p[idx];
   else if (idx%fN1==fN1-1) // C1 upper boundary
      E1[idx]=(V[idx-1]-V[idx])/dC1m[idx];
   else // bulk
      E1[idx]=(V[idx-1]-V[idx+1])/(dC1m[idx]+dC1p[idx]);

   return true;
}
//_____________________________________________________________________________
//
double X::GetC()
{
   Info("GetC","Start...");
   SuccessiveOverRelax(); // identify undepleted region
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
   SuccessiveOverRelax();
   // set impurity back
   if(fImpurity!=tmpImpurity) delete []fImpurity;
   fImpurity=tmpImpurity;

   // calculate C based on CV^2/2 = epsilon int E^2 dx^3 / 2
   double dV=Bias[0]-Bias[1]; if(dV<0)dV=-dV;
   double SumofElectricField=0;
   for(int i=0;i<fN;i++) {
      SumofElectricField+=E1[i]*E1[i]*dC1p[i]*cm*cm;
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
   bool b,d; double v,te,e1,e2,e3,c1,c2,c3;
   fTree = new TTree("t","field data");
   fTree->SetDirectory(0);
   fTree->Branch("v",&v,"v/D");
   fTree->Branch("e",&te,"e/D");
   // 1D data
   fTree->Branch("e1",&e1,"e1/D");
   fTree->Branch("c1",&c1,"c1/D");
   // initialize values
   if (dC1p[0]==0) Initialize(); // setup & initialize grid

   if (dC2p[0]!=0) { // if it is a 2D grid
      fTree->Branch("e2",&e2,"e2/D");
      fTree->Branch("c2",&c2,"c2/D");
   }
   if (dC3p[0]!=0) { // if it is a 3D grid
      fTree->Branch("e3",&e3,"e3/D");
      fTree->Branch("c3",&c3,"c3/D");
   }
   fTree->Branch("b",&b,"b/O"); // boundary flag
   fTree->Branch("d",&d,"d/O"); // depletion flag

   // fill tree
   Info("GetTree","%d entries",fN);
   for (int i=0; i<fN; i++) {
      e1= E1[i]; c1= C1[i]; // 1D data
      if (dC2p[i]!=0) { e2=E2[i]; c2=C2[i]; } // 2D data
      if (dC3p[i]!=0) { e3=E3[i]; c3=C3[i]; } // 3D data
      v = V[i]; b = fIsFixed[i]; d = fIsDepleted[i]; // common data
      if (dC3p[i]!=0) te=TMath::Sqrt(e1*e1 + e2*e2 + e3*e3);
      else { if (dC2p[i]!=0) te=TMath::Sqrt(e1*e1+e2*e2); else te=e1; }
      fTree->Fill();
   }

   fTree->GetListOfBranches()->ls();
   fTree->ResetBranchAddresses(); // disconnect from local variables
   return fTree;
}
//_____________________________________________________________________________
//
void X::SetGridImpurity()
{
   if (fImpDist && fImpurity[0]==0) // set impurity values if it's not done yet
      for (int i=fN;i-->0;) fImpurity[i]=fImpDist->Eval(C1[i], C2[i], C3[i]);
}
//_____________________________________________________________________________
//
int X::GetNsor()
{
   if (Gsor) return Gsor->GetX()[Gsor->GetN()-1];
   else return 0;
}

#include <TTree.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TStopwatch.h>

#include "Grid.h"
#include "Units.h"
using namespace GeFiCa;

FieldLine::FieldLine(const char *name, const char *title) : TNamed(name, title)
{
   // pick up a good style to modify
   gROOT->SetStyle("Plain");
   gStyle->SetName("GeFiCa");
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
Grid::Grid(size_t n1, size_t n2, size_t n3) : FieldLine("grid", "grid data"),
   N1(n1), N2(n2), N3(n3), RelaxationFactor(1.95), MaxIterations(5000),
   Precision(1e-7*volt), fTree(0)
{
}
//_____________________________________________________________________________
//
Grid::~X()
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
X& Grid::operator+=(GeFiCa::X *other)
{
   if (fN!=other->fN) {
      Warning("+=", 
            "Only same type of detector can be added together! Do nothing.");
      return *this; 
   }
   for (size_t i=0; i<fN; i++) {
      V[i]=V[i]+other->V[i];
      fImpurity[i]+=other->fImpurity[i];
   }
   Bias[0]+=other->Bias[0]; Bias[1]+=other->Bias[1]; 
   return *this;
}
//_____________________________________________________________________________
//
X& Grid::operator*=(double p)
{
   for (size_t i=0; i<fN; i++) V[i]=V[i]*p;
   Bias[0]*=p; Bias[1]*=p;
   return *this;
}
//_____________________________________________________________________________
//
size_t Grid::GetIdxOfMaxV()
{
   double max=V[0];
   size_t maxn=0;
   for(size_t i=1;i<fN;i++) {
      if(V[i]>max) {
         maxn=i;
         max=V[i];
      }
   }
   return maxn;
}
//_____________________________________________________________________________
//
size_t Grid::GetIdxOfMinV()
{
   double min=V[0];
   size_t minn=0;
   for(size_t i=1;i<fN;i++) {
      if(V[i]<min) {
         minn=i;
         min=V[i];
      }
   }
   return minn;
}
//_____________________________________________________________________________
//
bool Grid::IsDepleted()
{
   for(size_t i=0;i<fN;i++) {
      OverRelaxAt(i); // calculate one more time in case of 
      //adding two fields together, one is depleted, the other is not
      if (!fIsDepleted[i]) return false;
   }
   return true;
}
//_____________________________________________________________________________
//
void Grid::SetStepLength(double stepLength)
{
   for (size_t i=fN;i-->0;) {
      fIsFixed[i]=false;
      C1[i]=i*stepLength;
      dC1p[i]=stepLength;
      dC1m[i]=stepLength;
   }
}
//_____________________________________________________________________________
//
size_t* Grid::FindSurroundingMatrix(size_t idx)
{
   size_t *tmp=new size_t[3];
   tmp[0]=idx;
   if(idx-1<0)tmp[1]=1;
   else tmp[1]=idx-1;
   if(idx+1>=fN)tmp[2]=fN-2;
   else tmp[2]=idx+1;
   return tmp;
}
//_____________________________________________________________________________
//
bool Grid::SuccessiveOverRelax()
{
   if (dC1p[0]==0) Initialize(); // setup and initialize grid if it's not done

   Info("SuccessiveOverRelax","Start...");
   if (Gsor==0) {
      Gsor = new TGraph; Gsor->SetName("Gsor");
      Gsor->SetTitle(";Number of iterations;log10(precision)");
   }
   else Gsor->Set(0); // reset the graph
   double cp=1; // current presision
   size_t it=0; // # of iterations
   TStopwatch watch; watch.Start();
   while (it<MaxIterations) {
      if (it%100==0) {
         Printf("%4d steps, precision: %.1e (target: %.0e)", 
               it, cp, Precision);
         if (it!=0) Gsor->SetPoint(Gsor->GetN(),it,TMath::Log10(cp));
      }
      double XUpSum=0;
      double XDownSum=0;
      for (size_t i=0;i<fN;i++) {
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
   for (size_t i=0; i<fN; i++) if (!CalculateField(i)) return false;
   Printf("%4d steps, precision: %.1e (target: %.0e)", it, cp, Precision);
   Gsor->SetPoint(Gsor->GetN(),it,TMath::Log10(cp));
   Info("SuccessiveOverRelax", "CPU time: %.1f s", watch.CpuTime());
   return true;
}
//_____________________________________________________________________________
//
void Grid::OverRelaxAt(size_t idx)
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
size_t Grid::GetIdxOfPointNear(double c1,size_t begin,size_t end) const
{
   if (end==0) end=fN1-1;
   if (begin>=end)return end;
   size_t mid=(begin+end)/2;
   if(C1[mid]>=c1)return GetIdxOfPointNear(c1,begin,mid);
   else return GetIdxOfPointNear(c1,mid+1,end);
}
//_____________________________________________________________________________
//
size_t Grid::GetIdxOfPointNear(double c1,double c2,
      size_t begin,size_t end) const
{
   if (end==0) end=fN2-1;
   //search using binary search
   // if(begin>=end)cout<<"to x"<<begin<<" "<<end<<endl;;
   if(begin>=end)return GetIdxOfPointNear(c1,end*fN1,(end+1)*fN1-1);
   size_t mid=((begin+end)/2);
   if(C2[mid*fN1]>=c2){//cout<<"firsthalf"<<begin<<" "<<end<<endl; 
      return GetIdxOfPointNear(c1,c2,begin,mid);
   }
   else{//cout<<"senondhalf"<<begin<<" "<<end<<endl; 
      return GetIdxOfPointNear(c1,c2,mid+1,end);}
}
//_____________________________________________________________________________
//
size_t Grid::GetIdxOfPointNear(double c1, double c2,double c3,
      size_t begin,size_t end) const
{
   if (end==0) end=fN3-1;
   //search using binary search
   if(begin>=end)return GetIdxOfPointNear(c1,c2,begin,begin+fN1*fN2-1);
   size_t mid=((begin/(fN1*fN2)+end/(fN1*fN2))/2)*fN1*fN2;
   if(C3[mid]>=c3)return GetIdxOfPointNear(c1,c2,c3,begin,mid);
   else return GetIdxOfPointNear(c1,c2,c3,mid+1,end);
}
//_____________________________________________________________________________
//
double Grid::GetData(const vector<double> &data,
      double x, double y, double z) const
{
   size_t idx=GetIdxOfPointNear(x);
   if (idx==N1) return data[idx];
   double ab=(-x+C1[idx])/dC1p[idx];
   double aa=1-ab;
   return data[idx]*ab+data[idx-1]*aa;
}
//_____________________________________________________________________________
//
bool Grid::CalculateField(size_t idx)
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
double Grid::GetC()
{
   Info("GetC","Start...");
   SuccessiveOverRelax(); // identify undepleted region
   // set impurity to zero
   double *tmpImpurity=fImpurity;
   for (size_t i=0;i<fN;i++) {
      if (fImpurity[i]!=0) {
         fImpurity=new double[fN];
         for (size_t j=0;j<fN;j++) {
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
   for(size_t i=0;i<fN;i++) {
      SumofElectricField+=E1[i]*E1[i]*dC1p[i]*cm*cm;
      if (!fIsDepleted[i]) fIsFixed[i]=false;
   }
   double c=SumofElectricField*epsilon/dV/dV;
   Info("GetC","%.2f pF",c/pF);
   return c;
}
//_____________________________________________________________________________
//
TTree* Grid::GetTree(bool createNew)
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
   for (size_t i=0; i<fN; i++) {
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
void Grid::SetGridImpurity()
{
   if (fImpDist && fImpurity[0]==0) // set impurity values if it's not done yet
      for (size_t i=fN;i-->0;) fImpurity[i]=fImpDist->Eval(C1[i], C2[i], C3[i]);
}
//_____________________________________________________________________________
//
size_t Grid::GetNsor()
{
   if (Gsor) return Gsor->GetX()[Gsor->GetN()-1];
   else return 0;
}

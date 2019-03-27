#include <TTree.h>
#include <TGraph.h>

#include "Grid.h"
TGraph* FieldLine::GetGraph()
{
   if (GetN()<1) { Error("GetGraph", "No data for graph. Abort!"); abort(); }
   if (fGl) { return fGl; }
   fGl = new TGraph(GetN(),C1.data(),C2.data());
   fGl->SetName(GetName()); fGl->SetTitle(GetTitle());
   return fGl;
}
//_____________________________________________________________________________
//
#include <TStyle.h>
#include "Units.h"
using namespace GeFiCa;
Grid::Grid(size_t n1, size_t n2, size_t n3) : N1(n1), N2(n2), N3(n3),
   RelaxationFactor(1.95), MaxIterations(5000), Precision(1e-7*volt), fTree(0)
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
Grid& Grid::operator+=(Grid& other)
{
   if (GetN()!=other->GetN()) {
      Warning("+=", 
            "Numbers of points in two grids are different, return *this.");
      return *this; 
   }
   for (size_t i=0; i<GetN(); i++) {
      Vp[i]=Vp[i]+other.Vp[i];
      Src[i]=Src[i]+other->Src[i];
   }
   Bias[0]=Bias[0]+other->Bias[0]; Bias[1]=Bias[1]+other->Bias[1]; 
   return *this;
}
//_____________________________________________________________________________
//
Grid& Grid::operator*=(double scale)
{
   for (size_t i=0; i<GetN(); i++) Vp[i]=Vp[i]*scale;
   Bias[0]*=scale; Bias[1]*=scale;
   return *this;
}
//_____________________________________________________________________________
//
size_t Grid::GetIdxOfMaxV()
{
   double max=Vp[0]; size_t idx=0;
   for (size_t i=1; i<GetN(); i++) {
      if (Vp[i]>max) {
         idx=i;
         max=Vp[i];
      }
   }
   return idx;
}
//_____________________________________________________________________________
//
size_t Grid::GetIdxOfMinV()
{
   double min=Vp[0]; size_t idx=0;
   for (size_t i=1; i<GetN(); i++) {
      if (Vp[i]<min) {
         idx=i;
         min=Vp[i];
      }
   }
   return idx;
}
//_____________________________________________________________________________
//
bool Grid::IsDepleted()
{
   for(size_t i=0;i<GetN();i++) {
      // calculate one more time in case of adding two fields together,
      // one is depleted, the other is not
      OverRelaxAt(i);
      if (fIsDepleted[i]=false) return false;
   }
   return true;
}
//_____________________________________________________________________________
//
#include <TStopwatch.h>
void Grid::SuccessiveOverRelax()
{
   Info("SuccessiveOverRelax","Start...");
   double cp=1; // current presision
   size_t it=0; // # of iterations
   TStopwatch watch; watch.Start();
   while (it<MaxIterations) {
      if (it%100==0)
         Printf("%4d steps, precision: %.1e (target: %.0e)",it,cp,Precision);
      double numerator=0, denominator=0;
      for (size_t i=0; i<GetN(); i++) {
         double old=Vp[i]; // save old value of Vp[i]
         OverRelaxAt(i); // update Vp[i]
         if(old>0) denominator+=old; else denominator-=old;
         double diff=Vp[i]-old;
         if(diff>0) numerator+=(diff); else numerator-=(diff);
      }
      cp = numerator/denominator;
      it++;
      if (cp<Precision) break;
   }
   for (size_t i=0; i<GetN(); i++) CalculateE(i);
   Printf("%4d steps, precision: %.1e (target: %.0e)", it, cp, Precision);
   Info("SuccessiveOverRelax", "CPU time: %.1f s", watch.CpuTime());
}
//_____________________________________________________________________________
//
void Grid::OverRelaxAt(size_t idx)
{
   if (fIsFixed[idx])return ;
   double rho=-Src[idx]*Qe;
   double h2=dC1m[idx];
   double h3=dC1p[idx];
   double p2=Vp[idx-1];
   double p3=Vp[idx+1];

   double tmp=-rho/epsilon*h2*h3/2
      + (h3*Vp[idx-1]+h2*Vp[idx+1])/(h2+h3);

   //find min
   double min=p2;
   double max=p2;
   if(min>p3)min=p3;
   //find max
   if(max<p3)max=p3;
   //if tmp is greater or smaller than max and min, set tmp to it.

   //Vp[idx]=RelaxationFactor*(tmp-Vp[idx])+Vp[idx];
   double oldP=Vp[idx];
   tmp=RelaxationFactor*(tmp-oldP)+oldP;

   if(tmp<min) {
      Vp[idx]=min;
      fIsDepleted[idx]=false;
   } else if(tmp>max) {
      Vp[idx]=max;
      fIsDepleted[idx]=false;
   } else
      fIsDepleted[idx]=true;

   if(fIsDepleted[idx]||Bias[0]==Bias[1]) Vp[idx]=tmp;
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
void Grid::CalculateE(size_t idx)
{
   if (dC1p[idx]==0 || dC1m[idx]==0) return false;

   if (idx%fN1==0) // C1 lower boundary
      E1[idx]=(Vp[idx]-Vp[idx+1])/dC1p[idx];
   else if (idx%fN1==fN1-1) // C1 upper boundary
      E1[idx]=(Vp[idx-1]-Vp[idx])/dC1m[idx];
   else // bulk
      E1[idx]=(Vp[idx-1]-Vp[idx+1])/(dC1m[idx]+dC1p[idx]);
}
//_____________________________________________________________________________
//
double Grid::GetC()
{
   Info("GetC","Start...");
   SuccessiveOverRelax(); // identify undepleted region
   // set impurity to zero
   double *tmpImpurity=Src;
   for (size_t i=0;i<GetN();i++) {
      if (Src[i]!=0) {
         Src=new double[GetN()];
         for (size_t j=0;j<GetN();j++) {
            Src[j]=0;
            if (!fIsFixed[j] && !fIsDepleted[j]) fIsFixed[j]=true;
         }
         break;
      }
   }
   // calculate potential without impurity
   SuccessiveOverRelax();
   // set impurity back
   if(Src!=tmpImpurity) delete []Src;
   Src=tmpImpurity;

   // calculate C based on CV^2/2 = epsilon int E^2 dx^3 / 2
   double dV=Bias[0]-Bias[1]; if(dV<0)dV=-dV;
   double SumofElectricField=0;
   for(size_t i=0;i<GetN();i++) {
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
   Info("GetTree","%d entries",GetN());
   for (size_t i=0; i<GetN(); i++) {
      e1= E1[i]; c1= C1[i]; // 1D data
      if (dC2p[i]!=0) { e2=E2[i]; c2=C2[i]; } // 2D data
      if (dC3p[i]!=0) { e3=E3[i]; c3=C3[i]; } // 3D data
      v = Vp[i]; b = fIsFixed[i]; d = fIsDepleted[i]; // common data
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
   if (fImpDist && Src[0]==0) // set impurity values if it's not done yet
      for (size_t i=GetN();i-->0;) Src[i]=fImpDist->Eval(C1[i], C2[i], C3[i]);
}
//_____________________________________________________________________________
//
size_t Grid::GetNsor()
{
   if (Gsor) return Gsor->GetX()[Gsor->GetN()-1];
   else return 0;
}
//_____________________________________________________________________________
//

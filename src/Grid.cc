#include <TMath.h>
#include <TTree.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TStopwatch.h>

#include "Grid.h"
#include "Units.h"
#include "Detector.h"
using namespace GeFiCa;

TGraph* FieldLine::GetGraph()
{
   if (GetN()<1) {
      Warning("GetGraph", "No data for graph. Return 0!");
      return 0;
   }
   if (fGl) { return fGl; }
   fGl = new TGraph(GetN(),C1.data(),C2.data());
   fGl->SetName(GetName()); fGl->SetTitle(GetTitle());
   return fGl;
}
//______________________________________________________________________________
//

Grid::Grid(size_t n1, size_t n2, size_t n3) : N1(n1), N2(n2), N3(n3),
   MaxIterations(5000), RelaxationFactor(1.95), Precision(1e-7*volt),
   fTree(0), fDetector(0)
{
   // pick up a good style to modify
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
//______________________________________________________________________________
//
Grid& Grid::operator+=(Grid& other)
{
   if (GetN()!=other.GetN()) {
      Warning("+=", 
            "Numbers of points in two grids are different, return *this.");
      return *this; 
   }
   for (size_t i=0; i<GetN(); i++) {
      Vp[i]=Vp[i]+other.Vp[i];
      Src[i]=Src[i]+other.Src[i];
   }
   return *this;
}
//______________________________________________________________________________
//
Grid& Grid::operator*=(double scale)
{
   for (size_t i=0; i<GetN(); i++) Vp[i]=Vp[i]*scale;
   return *this;
}
//______________________________________________________________________________
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
//______________________________________________________________________________
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
//______________________________________________________________________________
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
//______________________________________________________________________________
//
void Grid::GetBoundaryConditionFrom(Detector &detector)
{
   if (GetN()>0) { // this function can only be called once
      Warning("GetBoundaryConditionFrom", "has been called. Do nothing.");
      return;
   }
}
//______________________________________________________________________________
//
void Grid::SuccessiveOverRelax()
{
   if (dC1p.size()<1) {
      Error("SuccessiveOverRelax", "Grid is not ready. "
            "Please call GetBoundaryConditionFrom(Detector&) first.");
      abort();
   }
   Info("SuccessiveOverRelax","Start...");
   double cp=1; // current presision
   fIterations=0;
   TStopwatch watch; watch.Start();
   while (fIterations<MaxIterations) {
      if (fIterations%100==0) Printf("%4zu steps, "
            "precision: %.1e (target: %.0e)", fIterations, cp, Precision);
      double numerator=0, denominator=0;
      for (size_t i=0; i<GetN(); i++) {
         double old=Vp[i]; // save old value of Vp[i]
         OverRelaxAt(i); // update Vp[i]
         if(old>0) denominator+=old; else denominator-=old;
         double diff=Vp[i]-old;
         if(diff>0) numerator+=(diff); else numerator-=(diff);
      }
      if (denominator==0) {
         Error("SuccessiveOverRelax","Sum of Vs == 0!"); abort();
      }
      cp = numerator/denominator;
      fIterations++;
      if (cp<Precision) break;
   }
   CalculateE();
   Printf("%4zu steps, precision: %.1e (target: %.0e)",
         fIterations, cp, Precision);
   Info("SuccessiveOverRelax", "CPU time: %.1f s", watch.CpuTime());
}
//______________________________________________________________________________
//
void Grid::SolveAnalytically()
{
   if (Src[0]!=Src[N1-1]) {
      Error("SolveAnalytically", "can't handle changing impurity.");
      abort();
   }
}
//______________________________________________________________________________
//
size_t Grid::GetIdxOfPointToTheRightOf(double c1,size_t begin,size_t end) const
{
   //if (end==0) end=N1-1;
   if (begin>=end)return end;
   size_t mid=(begin+end)/2;
   if(C1[mid]>=c1)return GetIdxOfPointToTheRightOf(c1,begin,mid);
   else return GetIdxOfPointToTheRightOf(c1,mid+1,end);
}
//______________________________________________________________________________
//
size_t Grid::GetIdxOfPointToTheRightOf(double c1,double c2,
      size_t begin,size_t end) const
{
   //if (end==0) end=N2-1;
   //search using binary search
   // if(begin>=end)cout<<"to x"<<begin<<" "<<end<<endl;;
   if(begin>=end)return GetIdxOfPointToTheRightOf(c1,end*N1,(end+1)*N1-1);
   size_t mid=((begin+end)/2);
   if(C2[mid*N1]>=c2) return GetIdxOfPointToTheRightOf(c1,c2,begin,mid);
   else return GetIdxOfPointToTheRightOf(c1,c2,mid+1,end);
}
//______________________________________________________________________________
//
size_t Grid::GetIdxOfPointToTheRightOf(double c1, double c2,double c3,
      size_t begin,size_t end) const
{
   //if (end==0) end=N3-1;
   //search using binary search
   if(begin>=end)return GetIdxOfPointToTheRightOf(c1,c2,begin,begin+N1*N2-1);
   size_t mid=((begin/(N1*N2)+end/(N1*N2))/2)*N1*N2;
   if(C3[mid]>=c3)return GetIdxOfPointToTheRightOf(c1,c2,c3,begin,mid);
   else return GetIdxOfPointToTheRightOf(c1,c2,c3,mid+1,end);
}
//______________________________________________________________________________
//
double Grid::GetC()
{
   Info("GetC","Start...");
   SuccessiveOverRelax(); // identify undepleted region
   std::vector<double>original=Src; // save original rho/epsilon
   for (size_t i=0; i<GetN(); i++) {
      Src[i]=0; // set impurity to zero
      if (fIsDepleted[i]==false) fIsFixed[i]=true; // fix undepleted points
   }
   SuccessiveOverRelax(); // calculate potential without impurity
   Src=original; // set impurity back

   return 0;
}
//______________________________________________________________________________
//
TTree* Grid::GetTree(bool createNew)
{
   if (fTree) { if (createNew) delete fTree; else return fTree; }

   // define tree
   bool b,d; double vp,et,e1,e2,e3,c1,c2,c3,p1,p2,p3,m1,m2,m3;
   fTree = new TTree("t","field data");
   fTree->SetDirectory(0);
   fTree->Branch("v",&vp,"v/D");
   fTree->Branch("e",&et,"e/D");
   // 1D data
   fTree->Branch("e1",&e1,"e1/D");
   fTree->Branch("c1",&c1,"c1/D");
   fTree->Branch("p1",&p1,"p1/D");
   fTree->Branch("m1",&m1,"m1/D");
   if (N2!=0) { // 2D data
      fTree->Branch("e2",&e2,"e2/D");
      fTree->Branch("c2",&c2,"c2/D");
      fTree->Branch("p2",&p2,"p2/D");
      fTree->Branch("m2",&m2,"m2/D");
   }
   if (N3!=0) { // 3D data
      fTree->Branch("e3",&e3,"e3/D");
      fTree->Branch("c3",&c3,"c3/D");
      fTree->Branch("p3",&p3,"p3/D");
      fTree->Branch("m3",&m3,"m3/D");
   }
   fTree->Branch("b",&b,"b/O"); // boundary flag
   fTree->Branch("d",&d,"d/O"); // depletion flag

   // fill tree
   Info("GetTree","%zu entries",GetN());
   for (size_t i=0; i<GetN(); i++) {
      e1=E1[i]; c1=C1[i]; p1=dC1p[i]; m1=dC1m[i]; // 1D data
      if (N2!=0) { e2=E2[i]; c2=C2[i]; p2=dC2p[i]; m2=dC2m[i]; } // 2D data
      if (N3!=0) { e3=E3[i]; c3=C3[i]; p3=dC3p[i]; m3=dC3m[i]; } // 3D data
      vp=Vp[i]; et=Et[i]; b=fIsFixed[i]; d=fIsDepleted[i]; // common data
      fTree->Fill();
   }

   fTree->GetListOfBranches()->ls();
   fTree->ResetBranchAddresses(); // disconnect from local variables
   return fTree;
}
//______________________________________________________________________________
//
double Grid::GetData(const std::vector<double> &data,
      double x, double y, double z) const
{
   if (N3!=0) { // 3D
      Warning("GetData", "not yet implemented. Abort."); abort();
   }
   if (N2!=0) { // 2D
      // https://codeplea.com/triangular-interpolation
      //       tmv
      // +-aa--+--ab--+(C1[idx], C2[idx])
      // |     ^      |
      // |  dyp|      ba
      // |     |(x,y) |
      // +<----+----->+
      // |     |      |
      // |  dym|      bb
      // |     v      |
      // +-----+------+
      //       bmv
      size_t idx=GetIdxOfPointToTheRightOf(x,y,0,N2-1); // always exists

      bool tl=false; // existence of top left grid point
      bool br=false; // existence of bottom right grid point
      bool bl=false; // existence of bottom left grid point
      if (idx%N1!=0) tl=true; // not left boundary
      if (idx>=N1) br=true; // not bottom boundary
      if (tl&&bl) bl=true; // neither left nor bottom boundary

      if(!tl&&!br&&!bl) { // bottom left corner
         return data[idx];
      } else if(tl&&!br&&!bl) { // bottom boundary
         return twopoint(new double[2] {data[idx-1],data[idx]},
               x,new double[2] {C1[idx-1],C1[idx]});
      } else if(!tl&&!bl&&br) { // left boundary
         return twopoint(new double[2]{data[idx-N1],data[idx]},
               y,new double[2]{C2[idx-N1],C2[idx]});
      }
      // no boundary case
      if(dC1p[idx-1]==dC1m[idx]&&dC1p[idx-N1-1]==dC1m[idx-N1]&&
            dC2p[idx-N1-1]==dC2m[idx-1]&&dC2p[idx-N1]==dC2m[idx]) {
         return fourpoint(new double[4]{data[idx-1],data[idx],data[idx-N1-1],data[idx-N1]},new double[2]{x,y},new double[4]{C1[idx-1],C1[idx],C1[idx-N1-1],C1[idx-N1]},new double[4]{C2[idx-1],C2[idx],C2[idx-N1-1],C2[idx-N1]});
      }

      //top left case o
      if(dC1p[idx-1]!=dC1m[idx]&&
            dC1p[idx-N1-1]==dC1m[idx-N1]&&
            dC2p[idx-N1-1]!=dC2m[idx-1]&&
            dC2p[idx-N1]==dC2m[idx])
      {
         double xb=dC1m[idx]<dC1p[idx-1] ? C1[idx]-dC1m[idx] : C1[idx-1]+dC1p[idx-1];
         double yb=dC2m[idx-1]<dC2p[idx-N1-1] ? C2[idx-1]-dC2m[idx-1] : C2[idx-N1-1]+dC2p[idx-N1-1];
         double bmv=(C1[idx-N1]-xb)/(C1[idx-N1]-C1[idx-N1-1])*data[idx-N1-1]+(xb-C1[idx-N1-1])/(C1[idx-N1]-C1[idx-N1-1])*data[idx-N1];
         double tmv=fIsFixed[idx]?data[idx]:data[idx-1];
         double mlv=fIsFixed[idx-1]?data[idx-1]:data[idx-N1-1];
         double mmv=(C2[idx]-yb)/(C2[idx]-C2[idx-N1])*bmv+(yb-C2[idx-N1])/(C2[idx]-C2[idx-N1])*tmv;
         if(x>xb)
         {
            return fourpoint(new double[4]{tmv,data[idx],bmv,data[idx-N1]},new double[2]{x,y},new double[4]{xb,C1[idx],xb,C1[idx-N1]},new double[4]{C2[idx],C2[idx],C2[idx-N1],C2[idx-N1]});
         }
         else if(y<yb)
         {
            return fourpoint(new double[4]{mlv,mmv,data[idx-N1-1],bmv},new double[2]{x,y},new double[4]{C1[idx-1],xb,C1[idx-1],xb},new double[4]{yb,yb,C2[idx-N1],C2[idx-N1]});
         }
         else if((C2[idx]-yb)/(xb-C1[idx-1]*(x-C1[idx-1])+yb>y))
         {
            return threepoint(new double[3]{tmv,mmv,mlv},new double[2]{x,y},new double[3]{xb,xb,C1[idx-1]},new double[3]{C2[idx],yb,yb});
         }
         else if((C2[idx]-yb)/(xb-C1[idx-1]*(x-C1[idx-1])+yb<y))
         {
            return threepoint(new double[3]{tmv,tmv,mlv},new double[2]{x,y},new double[3]{xb,C1[idx-1],C1[idx-1]},new double[3]{C2[idx],C2[idx],yb});
         }
      }
      //topright case o
      if(dC1p[idx-1]!=dC1m[idx]&&
            dC1p[idx-N1-1]==dC1m[idx-N1]&&
            dC2p[idx-N1-1]==dC2m[idx-1]&&
            dC2p[idx-N1]!=dC2m[idx])
      {
         double xb=dC1m[idx]<dC1p[idx-1] ? C1[idx]-dC1m[idx] : C1[idx-1]+dC1p[idx-1];
         double yb=dC2m[idx]<dC2p[idx-N1] ? C2[idx]-dC2m[idx] : C2[idx-N1]+dC2p[idx-N1];
         double bmv=(C1[idx-N1]-xb)/(C1[idx-N1]-C1[idx-N1-1])*data[idx-N1-1]+(xb-C1[idx-N1-1])/(C1[idx-N1]-C1[idx-N1-1])*data[idx-N1];
         double tmv=fIsFixed[idx]?data[idx]:data[idx-1];
         double mrv=fIsFixed[idx]?data[idx]:data[idx-N1];
         double mmv=(C2[idx]-yb)/(C2[idx]-C2[idx-N1])*bmv+(yb-C2[idx-N1])/(C2[idx]-C2[idx-N1])*tmv;
         if(x<xb)
         {
            return fourpoint(new double[4]{data[idx-1],tmv,data[idx-N1-1],bmv},new double[2]{x,y},new double[4]{C1[idx-1],xb,C1[idx-N1-1],xb},new double[4]{C2[idx],C2[idx],C2[idx-N1],C2[idx-N1]});
         }
         else if(y<yb)
         {
            return fourpoint(new double[4]{mmv,mrv,bmv,data[idx-N1]},new double[2]{x,y},new double[4]{xb,C1[idx],xb,C1[idx-1]},new double[4]{yb,yb,C2[idx-N1],C2[idx-N1]});
         }
         else if((C2[idx]-yb)/(C1[idx]-x*(x-C1[idx])+yb>y))
         {
            return threepoint(new double[3]{tmv,mmv,mrv},new double[2]{x,y},new double[3]{xb,xb,C1[idx]},new double[3]{C2[idx],yb,yb});
         }
         else if((C2[idx]-yb)/(C1[idx]-xb*(x-C1[idx])+yb<y))
         {
            return threepoint(new double[3]{tmv,data[idx],mrv},new double[2]{x,y},new double[3]{xb,C1[idx],C1[idx]},new double[3]{C2[idx],C2[idx],yb});
         }
      }
      //bottom left case o
      if(dC1p[idx-1]==dC1m[idx]&&
            dC1p[idx-N1-1]!=dC1m[idx-N1]&&
            dC2p[idx-N1-1]!=dC2m[idx-1]&&
            dC2p[idx-N1]==dC2m[idx])
      {
         double xb=dC1m[idx-N1]<dC1p[idx-N1-1] ? C1[idx-N1]-dC1m[idx-N1] : C1[idx-N1-1]+dC1p[idx-N1-1];
         double yb=dC2m[idx-1]<dC2p[idx-N1-1] ? C2[idx-1]-dC2m[idx-1] : C2[idx-N1-1]+dC2p[idx-N1-1];
         double tmv=(C1[idx]-xb)/(C1[idx]-C1[idx-1])*data[idx-1]+(xb-C1[idx-1])/(C1[idx]-C1[idx-1])*data[idx];
         double bmv=fIsFixed[idx]?data[idx]:data[idx-1];
         double mlv=fIsFixed[idx-1]?data[idx-1]:data[idx-N1-1];
         double mmv=(C2[idx]-yb)/(C2[idx]-C2[idx-N1])*bmv+(yb-C2[idx-N1])/(C2[idx]-C2[idx-N1])*tmv;
         if(x>xb)
         {
            return fourpoint(new double[4]{tmv,data[idx],bmv,data[idx-N1]},new double[2]{x,y},new double[4]{xb,C1[idx],xb,C1[idx-N1]},new double[4]{C2[idx],C2[idx],C2[idx-N1],C2[idx-N1]});
         }
         else if(y>yb)
         {
            return fourpoint(new double[4]{data[idx-N1-1],tmv,mlv,mmv},new double[2]{x,y},new double[4]{C1[idx-1],xb,C1[idx-1],xb},new double[4]{C2[idx],C2[idx],yb,yb});
         }
         else if((C2[idx]-yb)/(xb-C1[idx-1]*(x-C1[idx-1])+C2[idx-N1]>y))
         {
            return threepoint(new double[3]{bmv,mmv,mlv},new double[2]{x,y},new double[3]{xb,xb,C1[idx-1]},new double[3]{yb,yb,C2[idx-N1]});
         }
         else if((C2[idx]-yb)/(xb-C1[idx-1]*(x-C1[idx-1])+C2[idx-N1]<y))
         {
            return threepoint(new double[3]{data[idx-N1-1],bmv,mlv},new double[2]{x,y},new double[3]{C1[idx-1],xb,C1[idx-1]},new double[3]{C2[idx],C2[idx],yb});
         }
      }
      //bottom right case o
      if(dC1p[idx-1]==dC1m[idx]&&
            dC1p[idx-N1-1]!=dC1m[idx-N1]&&
            dC2p[idx-N1-1]==dC2m[idx-1]&&
            dC2p[idx-N1]!=dC2m[idx])
      {
         double xb=dC1m[idx-N1]<dC1p[idx-N1-1] ? C1[idx-N1]-dC1m[idx-N1] : C1[idx-N1-1]+dC1p[idx-N1-1];
         double yb=dC2m[idx]<dC2p[idx-N1] ? C2[idx]-dC2m[idx] : C2[idx-N1]+dC2p[idx-N1];
         double tmv=(C1[idx]-xb)/(C1[idx]-C1[idx-1])*data[idx-1]+(xb-C1[idx-1])/(C1[idx]-C1[idx-1])*data[idx];
         double bmv=fIsFixed[idx-N1]?data[idx-N1]:data[idx-N1-1];
         double mrv=fIsFixed[idx-N1]?data[idx-N1]:data[idx];
         double mmv=(C2[idx]-yb)/(C2[idx]-C2[idx-N1])*bmv+(yb-C2[idx-N1])/(C2[idx]-C2[idx-N1])*tmv;
         if(x<xb)
         {
            return fourpoint(new double[4]{data[idx-1],tmv,data[idx-N1-1],bmv},new double[2]{x,y},new double[4]{C1[idx-1],xb,C1[idx-N1-1],xb},new double[4]{C2[idx],C2[idx],C2[idx-N1],C2[idx-N1]});
         }
         else if(y>yb)
         {
            return fourpoint(new double[4]{tmv,data[idx],mmv,mrv},new double[2]{x,y},new double[4]{xb,C1[idx],xb,C1[idx]},new double[4]{C2[idx],C2[idx],yb,yb});
         }
         else if((yb-C2[idx-N1])/(C1[idx]-xb*(x-xb)+C2[idx-N1]>y))
         {
            return threepoint(new double[3]{bmv,mmv,mrv},new double[2]{x,y},new double[3]{xb,xb,C1[idx]},new double[3]{C2[idx-N1],yb,yb});
         }
         else if((yb-C2[idx-N1])/(C1[idx]-xb*(x-xb)+C2[idx-N1]<=y))
         {
            return threepoint(new double[3]{bmv,data[idx-N1],mrv},new double[2]{x,y},new double[4]{xb,C1[idx],C1[idx-1]},new double[4]{C2[idx-N1],C2[idx-N1],yb});
         }
      }
      //cross top down case o
      if(dC1p[idx-1]!=dC1m[idx]&&
            dC1p[idx-N1-1]!=dC1m[idx-N1]&&
            dC2p[idx-N1-1]==dC2m[idx-1]&&
            dC2p[idx-N1]==dC2m[idx])
      {
         double topb=dC1m[idx]<dC1p[idx-1] ? C1[idx]-dC1m[idx] : C1[idx-1]+dC1p[idx-1];
         double bottomb=dC1m[idx-N1]<dC1p[idx-N1-1] ? C1[idx-N1]-dC1m[idx-N1] : C1[idx-N1-1]+dC1p[idx-N1-1];
         double bmv=fIsFixed[idx-N1]?data[idx-N1]:data[idx-N1-1];
         double tmv=fIsFixed[idx]?data[idx]:data[idx-1];
         double TopValueWithBottomCorssPoint=topb>bottomb ? (topb-bottomb)/(topb-C1[idx-1])*data[idx-1]+(bottomb-C1[idx-1])/(topb-C1[idx-1])*tmv : 
            (bottomb-topb)/(C1[idx]-topb)*data[idx]+(C1[idx]-bottomb)/(C1[idx-1]-topb)*tmv ;
         double BottomValueWithTopCorssPoint=topb>bottomb ? (topb-bottomb)/(C1[idx-N1]-bottomb)*bmv+(C1[idx-N1]-topb)/(C1[idx-N1]-bottomb)*data[idx-N1] : 
            (topb-C1[idx-N1-1])/(bottomb-C1[idx-N1-1])*bmv+(bottomb-topb)/(bottomb-C1[idx-1-N1])*data[idx-N1-1] ;
         double TopLeft,TopRight,BottomLeft,BottomRight;
         double TopLeftV,TopRightV,BottomLeftV,BottomRightV;
         if(topb>bottomb) {
            TopLeft=bottomb;
            TopLeftV=TopValueWithBottomCorssPoint;
            TopRight=topb;
            TopRightV=tmv;
            BottomLeft=bottomb;
            BottomLeftV=bmv;
            BottomRight=topb;
            BottomRightV=BottomValueWithTopCorssPoint;
         } else {
            TopRight=bottomb;
            TopRightV=TopValueWithBottomCorssPoint;
            TopLeft=topb;
            TopLeftV=tmv;
            BottomRight=bottomb;
            BottomRightV=bmv;
            BottomLeft=topb;
            BottomLeftV=BottomValueWithTopCorssPoint;
         }
         if(x>TopRight) {
         if(x<1.185&&x>1.175&&y>4.915&&y<4.925)Printf("trv:%f,brv:%f",TopRightV,BottomRightV);
            return fourpoint(new double[4]{TopRightV,data[idx],BottomRightV,data[idx-N1]},
                  new double[2]{x,y},new double[4]{TopRight,C1[idx],BottomRight,C1[idx-N1]},
                  new double[4]{C2[idx],C2[idx],C2[idx-N1],C2[idx-N1]});
         } else if(x<TopLeft) {
            return fourpoint(new double[4]{data[idx-1],TopLeftV,data[idx-N1-1],BottomLeftV},
                  new double[2]{x,y},new double[4]{C1[idx-1],TopLeft,C1[idx-N1-1],BottomLeft},
                  new double[4]{C2[idx-1],C2[idx-1],C2[idx-N1-1],C2[idx-N1-1]});
         } else if((dC2m[idx])/((topb-bottomb)*(x-(TopRight-TopLeft)/2)+C2[idx-N1]>y-dC2m[idx]/2)) {
            return topb>bottomb ? threepoint(new double[3]{TopRightV,BottomLeftV,BottomRightV},new double[2]{x,y},new double[3]{TopRight,BottomLeft,BottomRight},new double[3]{C2[idx],C2[idx-N1],C2[idx-N1]}) :
               threepoint(new double[3]{TopLeft,BottomLeftV,BottomRightV},new double[2]{x,y},new double[3]{TopLeft,BottomLeft,BottomRight},new double[3]{C2[idx],C2[idx-N1],C2[idx-N1]});
         } else if((dC2m[idx])/((topb-bottomb)*(x-(TopRight-TopLeft)/2)+C2[idx-N1]<=y-dC2m[idx]/2)) {
            return topb<bottomb ? threepoint(new double[3]{TopRightV,TopLeftV,BottomRightV},new double[2]{x,y},new double[3]{TopRight,TopLeft,BottomRight},new double[3]{C2[idx],C2[idx],C2[idx-N1]}) :
               threepoint(new double[3]{TopLeft,TopLeftV,BottomLeftV},new double[2]{x,y},new double[3]{TopLeft,BottomLeft,BottomLeft},new double[3]{C2[idx],C2[idx],C2[idx-N1]});
         }
      }
      // cross left right case:
      //         +----+ ... C2[idx]
      //         |    |
      //     tlv + trv+ ... t2
      //         |\   |
      //         | \  |
      //         |  \ |
      //         |   \|
      //     blv + brv+ ... b2
      //         |    |
      //         |    |
      //         +----+ ... C2[idx-N1]
      //         .    .
      //         .    .
      // C1[idx-1]    C1[idx]
      if(dC1p[idx-1]==dC1m[idx]&&
            dC1p[idx-N1-1]==dC1m[idx-N1]&&
            dC2p[idx-N1-1]!=dC2m[idx-1]&&
            dC2p[idx-N1]!=dC2m[idx])
      {
         double leftb=dC2m[idx-1]<dC2p[idx-N1-1] ? C2[idx-1]-dC2m[idx-1] : C2[idx-N1-1]+dC2p[idx-N1-1]; // y of left intersection
         double rightb=dC2m[idx]<dC2p[idx-N1] ? C2[idx]-dC2m[idx] : C2[idx-N1]+dC2p[idx-N1]; // y of right intersection
         double lmv=fIsFixed[idx-1]?data[idx-1]:data[idx-N1-1]; // voltage of leftb
         double rmv=fIsFixed[idx]?data[idx]:data[idx-N1]; // voltage of rightb
         double LeftValueWithRightCorssPoint=leftb>rightb ?
            (leftb-rightb)/(leftb-C2[idx-N1-1])*data[idx-N1-1]+(rightb-C2[idx-N1-1])/(leftb-C2[idx-N1-1])*lmv : 
            (rightb-leftb)/(C2[idx-1]-leftb)*data[idx-1]+(C2[idx-1]-rightb)/(C2[idx-1]-leftb)*lmv;
         double RightValueWithLeftCorssPoint=leftb>rightb ?
            (rightb-leftb)/(C2[idx]-rightb)*data[idx]+(C2[idx]-leftb)/(C2[idx]-rightb)*rmv : 
            (leftb-C2[idx-N1])/(rightb-C2[idx-N1])*rmv+(rightb-leftb)/(rightb-C2[idx-N1])*data[idx-N1];
         double t2,b2,tlv,trv,blv,brv;

         if(leftb>rightb) {
            t2=leftb;
            tlv=lmv;
            trv=RightValueWithLeftCorssPoint;
            b2=rightb;
            blv=LeftValueWithRightCorssPoint;
            brv=rmv;
         } else {
            trv=rmv;
            t2=rightb;
            tlv=LeftValueWithRightCorssPoint;
            brv=RightValueWithLeftCorssPoint;
            b2=rightb;
            blv=lmv;
         }
         // find out whether the boundary goes from top-right to bottom-left or
         // from top-left to bottom-right and then do interpolation
         if(y>=t2) {
            return fourpoint(new double[4]{data[idx-1],data[idx],tlv,trv},
                  new double[2]{x,y},new double[4]{C1[idx-1],C1[idx],C1[idx-1],C1[idx]},
                  new double[4]{C2[idx],C2[idx],t2,t2});
         } else if(y<b2) {
            return fourpoint(new double[4]{blv,brv,data[idx-N1-1],data[idx-N1]},new double[2]{x,y},new double[4]{C1[idx-1],C1[idx],C1[idx-N1-1],C1[idx-N1]},new double[4]{b2,b2,C2[idx-N1-1],C2[idx-N1]});
         } else if((rightb-leftb)/(dC1m[idx])*(x-C1[idx-1])+(b2)>y) {
            return leftb>rightb ? threepoint(new double[3]{trv,blv,brv},new double[2]{x,y},new double[3]{C1[idx],C1[idx-N1],C1[idx-N1]},new double[3]{t2,b2,b2}) :
               threepoint(new double[3]{tlv,blv,brv},new double[2]{x,y},new double[3]{C1[idx],C1[idx-N1],C1[idx-N1]},new double[3]{t2,b2,b2});
         } else if((rightb-leftb)/(dC1m[idx])*(x-C1[idx-1])+(b2)<=y) {
            return leftb>rightb ? threepoint(new double[3]{tlv,blv,brv},new double[2]{x,y},new double[3]{C1[idx-1],C1[idx-1],C1[idx]},new double[3]{t2,b2,b2}) :
               threepoint(new double[3] {trv,blv,brv},new double[2]{x,y},new double[3]{C1[idx],C1[idx-1],C2[idx]},new double[3]{t2,b2,b2});
         }
      }
      return data[idx];
   } // 2D
   // 1D
   //     |<---dC1m[idx]--->|
   //     +---r1---+---r2---+
   // C1[idx-1]    x      C1[idx]
   size_t idx=GetIdxOfPointToTheRightOf(x,0,N1-1);
   double r2=(C1[idx]-x)/dC1m[idx];
   double r1=1-r2;
   double xval=data[idx]*r1+data[idx-1]*r2;
   if (gDebug>0) Info("GetData","data(x=%.4f)=%.2f, "
         "C1[%zu]=%.4f, C1[%zu]=%.4f, r1=%.2f, r2=%0.2f",
         x,xval,idx-1,C1[idx-1],idx,C1[idx],r1,r2);
   return xval;
}
//______________________________________________________________________________
//
void Grid::CalculateE()
{
   for (size_t i=0; i<GetN(); i++) { // deal with E1 only
      E1[i]=-(Vp[i+1]-Vp[i-1])/(dC1p[i]+dC1m[i]); Et[i]=E1[i];
   }
   E1[0]=-(Vp[1]-Vp[0])/dC1p[0]; Et[0]=E1[0];
   E1[N1-1]=-(Vp[N1-1]-Vp[N1-2])/dC1m[N1-1]; Et[N1-1]=E1[N1-1];
}
//______________________________________________________________________________
//
double Grid::twopoint(double dataset[2],double tarlocationset,
      double pointxset[2])const
{
   double ab=(tarlocationset-pointxset[0])/(pointxset[1]-pointxset[0]);

   double result=dataset[0]*(1-ab) + dataset[1]*ab;
   delete [] dataset;
   delete [] pointxset;
   if (TMath::Abs(result)<1e-10) return 0;
   return result;
}
//______________________________________________________________________________
//
double Grid::threepoint(double dataset[3],double tarlocationset[2],
      double pointxset[3],double pointyset[3])const
{
   double x=tarlocationset[0];
   double x1=pointxset[0];
   double x2=pointxset[1];
   double x3=pointxset[2];

   double y=tarlocationset[1];
   double y1=pointyset[0];
   double y2=pointyset[1];
   double y3=pointyset[2];
   //x=(1-u-v)*x1+u*x2+v*x3
   //y=(1-u-v)*y1+u*y2+v*y3
   double v=((y2-y1)*(x-x1)-(y-y1)*(x2-x1))/((x3-x1*(y2-y1)-(y3-y1)*(x2-x1)));
   double u=((y3-y1)*(x-x1)-(y-y1)*(x3-x1))/((x2-x1*(y3-y1)-(y2-y1)*(x3-x1)));
   //double u=(y-y1-v*(y3-y1))/(y2-y1);

   double result=dataset[0]*(1-u-v)+dataset[1]*u+dataset[2]*v;

   delete [] dataset;
   delete [] tarlocationset;
   delete [] pointxset;
   delete [] pointyset;

   if (TMath::Abs(result)<1e-10) return 0;
   return result;
}
//______________________________________________________________________________
//
double Grid::fourpoint(double dataset[4],double tarlocationset[2],
      double pointxset[4],double pointyset[4])const
{
   //0---------1
   //|         |
   //|         |
   //|         |
   //2---------3
   double ab=(pointxset[1]-tarlocationset[0])/(pointxset[1]-pointxset[0]);
   double aa=1-ab;
   double ba=(tarlocationset[1]-pointyset[3])/(pointyset[1]-pointyset[3]);
   double bb=1-ba;

   double result=(dataset[1]*aa+dataset[0]*ab)*ba+(dataset[2]*ab+dataset[3]*aa)*bb;
   delete [] dataset;
   delete [] tarlocationset;
   delete [] pointxset;
   delete [] pointyset;

   if (TMath::Abs(result)<1e-10) return 0;
   return result;
}

//______________________________________________________________________________
//
FieldLine* Grid::GetFieldLineFrom(double x, double y, double z, bool positive)
{
   FieldLine *line=new FieldLine();

   int i=0;
   while (true) { // if (x,y) is in crystal
      line->C1.push_back(x);
      line->C2.push_back(y);
      if (x>=C1[GetN()-1]||x<=C1[0]||y>=C2[GetN()-1]||y<=C2[0]
            || fIsFixed[GetIdxOfPointToTheRightOf(x,y,0,N2)]) {//out of crystal
         if (i==0) // initial point is not in crystal
            Warning("GetFieldLineFrom", "Start point (%.1fcm,%.1fcm)"
                  " is not in crystal! Stop propagating.", x/cm, y/cm);
         break;
      }

      double ex = GetE1(x,y), ey = GetE2(x,y);
      double et = TMath::Sqrt(ex*ex+ey*ey); // total E
      if (et==0) {
         Warning("GetFieldLineFrom", "E@(x=%.1fmm,y=%.1fmm)=%.1fV/cm!",
               x/mm, y/mm, et/volt*cm);
         break;
      }
      double weight=5/(et/volt*mm); // propagate more in weaker field
      double dt=10*ns, mu=50000*cm2/volt/sec; // mu is mobility
      double dx=mu*ex*dt*weight, dy=mu*ey*dt*weight;
      line->dC1p.push_back(dx);
      line->dC2p.push_back(dy);
      if (gDebug>0)
         Printf("%04d x=%.3fmm, y=%.3fmm, Ex=%.2fV/cm, Ey=%.2fV/cm, "
               "weight=%.3f, dx=%.3fmm, dy=%.3fmm", i, x/mm, y/mm, ex/volt*cm,
               ey/volt*cm, weight, dx/mm, dy/mm);
      if (i>2000) {
         Info("GetFieldLineFrom", "Propagated more than 2000 steps. Stop");
         break;
      }
      if (positive) { x+=dx; y+=dy; } else { x-=dx; y-=dy; }
      i++;
   }

   return line;
}

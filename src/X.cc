#include <TFile.h>
#include <TChain.h>
#include <TVectorD.h>
#include <TStopwatch.h>
#include <TF3.h>
#include <Math/Functor.h>

#include "X.h"
#include "Units.h"
using namespace GeFiCa;
using namespace std;

X::X(int nx, const char *name, const char *title) : TNamed(name,title), n1(nx),
   n(nx), Csor(1.95), Precision(1e-7), MaxIterations(1e5), V0(0), V1(2e3*volt)
{
   if (n<10) { Warning("X","n<10, set it to 11"); n=11; n1=11; }

   fV=new double[n];
   fE1=new double[n];
   fC1=new double[n];
   fdC1p=new double[n];
   fdC1m=new double[n];
   fIsFixed=new bool[n];
   fIsDepleted=new bool[n];
   fImpurity=new double[n];

   for (int i=0;i<n;i++) {
      fV[i]=0;
      fE1[i]=0;
      fC1[i]=0;
      fdC1m[i]=0;
      fdC1p[i]=0;
      fIsFixed[i]=false;
      fIsDepleted[i]=true;
      fImpurity[i]=0;
   }

   fTree=NULL;
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
   if (method==kAnalytic) return Analytic();

   cout<<" Calculate field ..."<<endl;
   TStopwatch watch; watch.Start();
   int cnt=0;
   while (cnt++<MaxIterations) {
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
      double cp = XUpSum/XDownSum; // current precision
      if (cnt%100==0)
         Printf("  %05d iterations, precision: %e (target: %.0e)", 
               cnt, cp, Precision);

      if (cp<Precision) {
         for (int i=0; i<n; i++) if (!CalculateField(i)) return false;
         Printf("  %05d iterations, precision: %e (target: %.0e)", 
               cnt, cp, Precision);
         cout<<" Done. Spent "; watch.Stop(); watch.Print();
         return true;
      }
   }
   return false;
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
void X::SaveField(const char *fout)
{
   Info("SaveField", "%s", fout);

   TFile *file = new TFile(fout,"recreate");
   TTree *tree = new TTree("t","field tree"); // create it inside the file

   // variables
   TVectorD var(10);
   var[0]=(double)MaxIterations;
   var[1]=(double)n; // n=n1xn2xn3
   var[2]=Csor;
   var[7]=(double)n1;
   var[8]=1; // n2
   var[9]=1; // n3
   var.Write("v");

   // fields
   bool fix,dep; double e1,c1,v,dc1p,dc1m,im;
   tree->Branch("e1",&e1,"e1/D"); // electric field in the 1st coordinate
   tree->Branch("c1",&c1,"c1/D"); // 1st coordinate
   tree->Branch("v",&v,"v/D"); // electric potential
   tree->Branch("dc1p",&dc1p,"dc1p/D"); // step length to next point in c1
   tree->Branch("dc1m",&dc1m,"dc1m/D"); // step length to previous point in c1
   tree->Branch("fix",&fix,"fix/O"); // are the values fixed on this point
   tree->Branch("dep",&dep,"dep/O"); // is this point depleted
   tree->Branch("im",&im,"im/D"); // impurity

   for(int i=0;i<n;i++) {
      v = fV[i];
      e1 =fE1[i];
      c1 =fC1[i];
      im =fImpurity[i];
      fix=fIsFixed[i];
      dep=fIsDepleted[i];
      dc1p=fdC1p[i];
      dc1m=fdC1m[i];

      tree->Fill();
   }
   file->Write();
   file->Close(); // tree is deleted here
   delete file;
}
//_____________________________________________________________________________
//
void X::LoadField(const char *fin)
{
   Info("LoadField", "%s", fin);
   TFile *file = new TFile(fin);
   if (file->IsZombie()) Fatal("LoadField", "%s cannot be loaded", fin);

   // load list of variables
   TVectorD *variables = (TVectorD*) file->Get("v");
   double *var = variables->GetMatrixArray();
   MaxIterations = (int)var[0];
   n = (int)var[1];
   Csor = var[2];
   n1	= (int)var[7];

   file->Close();
   delete file;

   // load fields
   TChain *t = new TChain("t");
   t->Add(fin);
   bool fix,dep; double e1,c1,v,dc1p,dc1m,im;
   t->SetBranchAddress("e1",&e1);
   t->SetBranchAddress("c1",&c1);
   t->SetBranchAddress("v",&v);
   t->SetBranchAddress("dc1p",&dc1p);
   t->SetBranchAddress("dc1m",&dc1m);
   t->SetBranchAddress("im",&im);
   t->SetBranchAddress("fix",&fix);
   t->SetBranchAddress("dep",&dep);

   this->~X(); // delete old arrays if there are
   fV = new double[n];
   fE1 = new double[n];
   fC1 = new double[n];
   fdC1p = new double[n];
   fdC1m = new double[n];
   fImpurity = new double[n];
   fIsFixed = new bool[n];
   fIsDepleted = new bool[n];

   for (int i=0;i<n;i++) {
      t->GetEntry(i);
      fV[i]=v;
      fE1[i]=e1;
      fC1[i]=c1;
      fdC1p[i]=dc1p;
      fdC1m[i]=dc1m;
      fImpurity[i]=im;
      fIsFixed[i]=fix;
      fIsDepleted[i]=dep;
   }

   delete t;
}
//_____________________________________________________________________________
//
void X::SetImpurity(TF3 *fi)
{
   Initialize();
   for (int i=n;i-->0;) fImpurity[i]=fi->Eval(fC1[i]);
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
//_____________________________________________________________________________
//
TTree* X::GetTree()
{
   if (fTree!=NULL) return fTree;

   bool b,d; double e1,c1,v,dc1p,dc1m,i;
   // define tree
   fTree=new TTree("t","field data");
   fTree->Branch("v",&v,"v/D");
   fTree->Branch("e1",&e1,"e1/D");
   fTree->Branch("c1",&c1,"c1/D");
   fTree->Branch("dc1p",&dc1p,"dc1p/D"); // step length to next point
   fTree->Branch("dc1m",&dc1m,"dc1m/D"); // step length to previous point
   fTree->Branch("b",&b,"b/O"); // boundary flag
   fTree->Branch("d",&d,"d/O"); // depletion flag
   fTree->Branch("i",&i,"i/D"); // impurity
   // fill tree
   if (fdC1p[0]==0) Initialize(); // setup & initialize grid
   for(int j=0;j<n;j++) {
      v = fV[j];
      e1= fE1[j];
      c1= fC1[j];
      i = fImpurity[j];
      b = fIsFixed[j];
      d = fIsDepleted[j];
      dc1p=fdC1p[j];
      dc1m=fdC1m[j];

      fTree->Fill();
   }

   fTree->ResetBranchAddresses(); // disconnect from local variables
   return fTree;
}

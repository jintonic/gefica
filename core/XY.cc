#include <TF2.h>
#include <TFile.h>
#include <TChain.h>
#include <TVectorD.h>

#include "XY.h"
using namespace GEFICA;
//______________________________________________________________________________
// Create grid for 2-D field calculation.
ClassImp(XY)
XY::XY(unsigned short nx, unsigned short ny): X(nx*ny), n2(ny),
   fE2(0), fC2(0), fDistanceToLeft(0), fDistanceToRight(0)
{
  //claim a 2D field with n1*n2 Grid
   n=nx*ny; 
   n1=nx;
   fE2=new double[n];
   fC2=new double[n];
   fDistanceToLeft=new double[n];
   fDistanceToRight=new double[n];}

XY::~XY()
{
   if (fE2) delete[] fE2;
   if (fC2) delete[] fC2;
   if (fDistanceToLeft) delete[] fDistanceToLeft;
   if (fDistanceToRight) delete[] fDistanceToRight;
}

void XY::SetStepLength(double steplength1,double steplength2)
{
//set field step length
   X::SetStepLength(steplength1);
   for (int i=0;i<n;i++) {
      if(i>n1-1)fC2[i]=fC2[i-n1]+steplength2;
      else fC2[i]=0;
      if(i%n1==0)fC1[i]=0;
      else fC1[i]=fC1[i-1]+steplength1;

      fE2[i]=0;
      fDistanceToLeft[i]=steplength2;
      fDistanceToRight[i]=steplength2;
   }
}

void XY::SOR2(int idx,bool elec)
{
  
  // 2nd-order Runge-Kutta Successive Over-Relaxation
   if (fIsFixed[idx])return;
   double density=fImpurity[idx]*1.6e12;
   double h2=fDistanceToPrevious[idx];
   double h3=fDistanceToNext[idx];
   double h4=fDistanceToLeft[idx];
   double h1=fDistanceToRight[idx];
   double Pym1,Pyp1,Pxm1,Pxp1;
   if(idx>=n1)Pym1=fPotential[idx-n1];
   else Pym1=fPotential[idx];
   if(idx>=n-n1)Pyp1=fPotential[idx];
   else Pyp1=fPotential[idx+n1];
   if(idx%n1==0)Pxm1=fPotential[idx];
   else Pxm1=fPotential[idx-1];
   if(idx%n1==n1-1)Pxp1=fPotential[idx];
   else Pxp1=fPotential[idx+1];
   double tmp=((h1+h4)*(h1*h2*h4*Pxp1+h1*h3*h4*Pxm1)+(h2+h3)*(h1*h2*h3*Pyp1+h2*h3*h4*Pym1)-0.5*density/epsilon*(h1+h4)*(h2+h3)*h1*h2*h3*h4)/((h1+h4)*(h1*h2*h4+h1*h3*h4)+(h2+h3)*(h1*h2*h3+h2*h3*h4));
   fPotential[idx]=Csor*(tmp-fPotential[idx])+fPotential[idx];
   if(elec)
   {
     fE1[idx]=(Pxp1-Pxm1)/(h2+h3);
     fE2[idx]=(Pyp1-Pym1)/(h1+h4);
   }
}

int XY::FindIdx(double tarx,double tary ,int ybegin,int yend)
{
  //search using binary search
   if(ybegin>=yend)return X::FindIdx(tarx,ybegin,ybegin+n1-1);
   int mid=((ybegin/n1+yend/n1)/2)*n1;
   if(fC2[mid]>=tary)return FindIdx(tarx,tary,ybegin,mid);
   else return FindIdx(tarx,tary,mid+1,yend);
}

double XY::GetData(double tarx, double tary, int thing)
{
  //ask thing with coordinate and item number: 1:Impurity 2:Potential 3:E1 4:E2

   int idx=FindIdx(tarx,tary,0,n);
   double ab=(tarx-fC1[idx])/fDistanceToNext[idx];
   double aa=1-ab;
   double ba=(tary-fC2[idx])/fDistanceToRight[idx];
   double bb=1-ba;
   double tar0,tar1,tar2,tar3,*tar=NULL;
   switch(thing)
   {
      case 0:tar= fImpurity;break;
      case 1:tar= fPotential;break;
      case 2:tar= fE1;break;
      case 3:tar= fE2;break;
   }
   tar3=-1;
   tar0=tar[idx];
   if(idx/n1+1==n1){tar1=0;tar3=0;}
   else {tar1=tar[idx+1];}
   if(idx>n-n1){tar2=0;tar3=0;}
   else {tar2=tar[idx+n1];}
   if (tar3==-1)tar3=tar[idx+n1+1];
   return (tar0*aa+tar1*ab)*ba+(tar2*aa+tar3*ab)*bb;
}
void XY::SaveField(const char * fout)
{
   X::SaveField(fout);
   TFile *file=new TFile(fout,"update");
   TVectorD  v=*(TVectorD*)file->Get("v");
   v[8]=(double)n2;
   v.Write();
   TTree * tree=(TTree*)file->Get("t");
   double E2s,C2s;
   tree->Branch("e2",&E2s,"e2/D"); // Electric field in y
   tree->Branch("c2",&C2s,"c2/D"); // persition in y
   for(int i=0;i<n;i++) {
      E2s=fE2[i];
      C2s=fC2[i];
      tree->Fill();
   }
   file->Write();
   file->Close();
   delete file;

}
void XY::LoadField(const char * fin)
{
  //will calculate electric field after load
   X::LoadField(fin);
   TFile *file=new TFile(fin);
   TVectorD *v1=(TVectorD*)file->Get("v");
   double * v=v1->GetMatrixArray();
   n2		=(int)	v[8];

   TChain *t =new TChain("t");
   t->Add(fin);
   double fEy,fPy;
   t->SetBranchAddress("c2",&fPy);
   t->SetBranchAddress("e2",&fEy);

   fE2=new double[n];
   fC2=new double[n];

   for (int i=0;i<n;i++) {
      t->GetEntry(i);
      fE2[i]=fEy;
      fC2[i]=fPy;
   }
   file->Close();
   delete file;
}
void XY::SetImpurity(TF2 * Im)
{
   for(int i=n;i-->0;) {
      fImpurity[i]=Im->Eval((double)fC1[i],(double)fC2[i]);
   }
}
void XY::SetImpurity(double density)
{
   for(int i=n;i-->0;) fImpurity[i]=density;
}

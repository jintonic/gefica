#include <TF3.h>
#include <TFile.h>
#include <TChain.h>
#include <TVectorD.h>

#include "XYZ.h"
using namespace GEFICA;

XYZ::XYZ(unsigned short n1, unsigned short n2,unsigned short n3): 
   XY(n1,n2), n3(n3), fE3(0), fC3(0), fDistanceToUp(0), fDistanceToDown(0)
{n=n1*n2*n3;}

XYZ::~XYZ()
{
   if (fE3) delete[] fE3;
   if (fC3) delete[] fC3;
   if (fDistanceToUp) delete[] fDistanceToUp;
   if (fDistanceToDown) delete[] fDistanceToDown; 
}

void XYZ::CreateGridWithFixedStepLength(double steplength)
{
   XY::CreateGridWithFixedStepLength(steplength);
   fE3=new double[n];
   fC3=new double[n];
   fDistanceToUp=new double[n];
   fDistanceToDown=new double[n];
   for (int i=0;i<n;i++) {
      if(i/n1*n2==0)fC3[i]=0;
      else fC3[i]=fC3[i-9]+steplength;
      if((i%(n1*n2))/n1!=0)fC2[i]=fC2[i-n1]+steplength;
      else fC2[i]=0;
      if(i%n1==0)fC1[i]=0;
      else fC1[i]=fC1[i-1]+steplength;

      fE3[i]=0;
      fDistanceToLeft[i]=steplength;
      fDistanceToRight[i]=steplength;
   }
}

void XYZ::Update(int idx)
{//need update
   if (fIsFixed[idx])return;
   double density=fImpurity[idx]*1.6e12;
   double h2=fDistanceToPrevious[idx];
   double h3=fDistanceToNext[idx];
   double h4=fDistanceToLeft[idx];
   double h1=fDistanceToRight[idx];
   double h0=fDistanceToDown[idx];
   double h5=fDistanceToUp[idx];
   double Pym1,Pyp1,Pxm1,Pxp1,Pzp1,Pzm1;
   if(idx<n1*n2)Pzm1=fPotential[idx];
   else Pzm1=fPotential[idx-n1*n2];
   if(idx>=n-n1*n2)Pzp1=fPotential[idx];
   else Pzp1=fPotential[idx+n1*n2];
   if(idx%(n1*n2)>(n1*n2)-n1-1) Pyp1=fPotential[idx];
   else Pyp1=fPotential[idx+n1];
   if(idx%(n1*n2)<n1)Pym1=fPotential[idx];
   else Pym1=fPotential[idx-n1];
   if((idx%(n1*n2))%n1==n1-1)Pxp1=fPotential[idx];
   else Pxp1=fPotential[idx+1];
   if((idx%(n1*n2))%n1==0)Pxm1=fPotential[idx];
   else Pxm1=fPotential[idx-1];

   double tmp= (
         -density/epsilon*h0*h1*h2*h3*h4*h5*(h1+h4)*(h2+h3)*(h0+h5)
         +(Pxp1*h3+Pxm1*h2)*h0*h1*h4*h5*(h1+h4)*(h0+h5)
         +(Pyp1*h4+Pym1*h1)*h0*h2*h3*h5*(h0+h5)*(h2+h3)
         +(Pzp1*h5+Pzm1*h0)*h1*h2*h3*h4*(h1+h4)*(h2+h3)	
         )
      /((h0+h5)*(h1+h4)*(h2+h3)*(h0*h1*h4*h5+h0*h2*h3*h5+h1*h2*h3*h4));

   fPotential[idx]=Csor*(tmp-fPotential[idx])+fPotential[idx];
   fE1[idx]=(Pxp1-Pxm1)/(h2+h3);
   fE2[idx]=(Pyp1-Pym1)/(h1+h4);
   fE3[idx]=(Pzp1-Pzm1)/(h0+h5);
}

int XYZ::FindIdx(double tarx, double tary ,double tarz,int begin,int end)
{
   if(begin>=end)return XY::FindIdx(tarx,tary,begin,begin+n1*n2-1);
   int mid=((begin/(n1*n2)+end/(n1*n2))/2)*n1*n2;
   if(fC3[mid]>=tarz)return FindIdx(tarx,tary,tarz,begin,mid);
   else return FindIdx(tarx,tary,tarz,mid+1,end);
}

double XYZ::GetData(double tarx, double tary, double tarz,int thing)
{
   int idx=FindIdx(tarx,tary,tarz,0,n);
   double ab=(tarx-fC1[idx])/fDistanceToNext[idx];
   double aa=1-ab;
   double ba=(tary-fC2[idx])/fDistanceToRight[idx];
   double bb=1-ba;
   double ac=(tarz-fC3[idx])/fDistanceToUp[idx];
   double ca=1-ac;
   double tar0,tar1,tar2,tar3,tar4,tar5,tar6,tar7,*tar=NULL;
   switch(thing)
   {
      case 0:tar= fImpurity;break;
      case 1:tar= fPotential;break;
      case 2:tar= fE1;break;
      case 3:tar= fE2;break;
      case 4:tar= fE3;break;
   }
   tar3=-1;
   tar5=-1;
   tar6=-1;
   tar7=-1;
   tar0=tar[idx];
   if(idx>=(n-n1*n2)){tar4=0;tar5=0;tar6=0;tar7=0;}
   else{tar4=tar[idx+n1*n2];}
   if(idx%(n1*n2)%n1==n1-1){tar2=0;tar3=0;tar6=0;tar7=0;}
   else{tar2=tar[idx+n1];}
   if(idx%(n1*n2)/n1==n2-1){tar1=0;tar3=0;tar5=0;tar7=0;}
   if(tar3==-1)tar3=tar[idx+n1+1];
   if(tar5==-1)tar5=tar[idx+n1*n2+1];
   if(tar6==-1)tar6=tar[idx+n1*n2+n1];
   if(tar7==-1)tar7=tar[idx+n1*n2+n1+1];
   return ((tar0*aa+tar1*ab)*ba+(tar2*aa+tar3*ab)*bb)*ac+((tar4*aa+tar5*ab)*ba+(tar6*aa+tar7*ab)*bb)*ca;
}
void XYZ::Save(const char * fout)
{
   XY::Save(fout);
   TFile *file=new TFile(fout,"update");
   TVectorD  v=*(TVectorD*)file->Get("v");
   v[9]=(double)n3;
   v.Write();
   TTree * tree=(TTree*)file->Get("t");
   double E3s,C3s;
   tree->Branch("e3",&E3s,"e3/D"); // Electric field in z
   tree->Branch("c3",&C3s,"c3/D"); // persition in z
   for(int i=0;i<n;i++) {
      E3s=fE3[i];
      C3s=fC3[i];
      tree->Fill();
   }
   file->Write();
   file->Close();
   delete file;

}
void XYZ::Load(const char * fin)
{
   XY::Load(fin);
   TFile *file=new TFile(fin);
   TVectorD *v1=(TVectorD*)file->Get("v");
   double * v=v1->GetMatrixArray();
   n3 = (int) v[9];

   TChain *t =new TChain("t");
   t->Add(fin);
   double fEz,fPz;
   t->SetBranchAddress("c3",&fPz);
   t->SetBranchAddress("e3",&fEz);

   fE3=new double[n];
   fC3=new double[n];

   for (int i=0;i<n;i++) {
      t->GetEntry(i);
      fE3[i]=fEz;
      fC3[i]=fPz;
   }
   file->Close();
   delete file;
}
void XYZ::SetImpurity(TF3 * Im)
{
   for(int i=n;i-->0;) {
      fImpurity[i]=Im->Eval((double)fC1[i],(double)fC2[i],(double)fC3[i]);
   }
}

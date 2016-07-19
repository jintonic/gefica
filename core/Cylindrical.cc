#include <TF3.h>
#include <TFile.h>
#include <TChain.h>
#include <TVectorD.h>

#include "Cylindrical.h"
using namespace GEFICA;

Cylindrical::~Cylindrical()
{
   delete[] E1;
   delete[] E2;
   delete[] E3;
   delete[] C1;
   delete[] C2;
   delete[] C3;
   delete[] P;
   delete[] StepNext;
   delete[] StepBefore;
   delete[] StepLeft;
   delete[] StepRight;
   delete[] StepUp;
   delete[] StepDown;
   delete[] isbegin;
   delete[] Impurity;
}

void Cylindrical::Create(double steplength)
{
   XY::Create(steplength);
   E3=new double[n];
   C3=new double[n];
   StepUp=new double[n];
   StepDown=new double[n];
   for (int i=0;i<n;i++) {
      if(i/x*y==0)C3[i]=0;
      else C3[i]=C3[i-9]+steplength;
      if((i%(x*y))/x!=0)C2[i]=C2[i-x]+steplength;
      else C2[i]=0;
      if(i%x==0)C1[i]=0;
      else C1[i]=C1[i-1]+steplength;

      E3[i]=0;
      StepLeft[i]=steplength;
      StepRight[i]=steplength;
   }
}

void Cylindrical::Update(int idx)
{//need update
   if (isbegin[idx])return;
   double density=Impurity[idx]*1.6e12;
   double h2=StepBefore[idx];
   double h3=StepNext[idx];
   double h4=StepLeft[idx];
   double h1=StepRight[idx];
   double h0=StepDown[idx];
   double h5=StepUp[idx];
   double Pym1,Pyp1,Pxm1,Pxp1,Pzp1,Pzm1;
   if(idx<x*y)Pzm1=P[idx];
   else Pzm1=P[idx-x*y];
   if(idx>=n-x*y)Pzp1=P[idx];
   else Pzp1=P[idx+x*y];
   if(idx%(x*y)>(x*y)-x-1) Pyp1=P[idx];
   else Pyp1=P[idx+x];
   if(idx%(x*y)<x)Pym1=P[idx];
   else Pym1=P[idx-x];
   if((idx%(x*y))%x==x-1)Pxp1=P[idx];
   else Pxp1=P[idx+1];
   if((idx%(x*y))%x==0)Pxm1=P[idx];
   else Pxm1=P[idx-1];

   double tmp= (
         -density/(E0*ER)*h0*h1*h2*h3*h4*h5*(h1+h4)*(h2+h3)*(h0+h5)
         +(Pxp1*h3+Pxm1*h2)*h0*h1*h4*h5*(h1+h4)*(h0+h5)
         +(Pyp1*h4+Pym1*h1)*h0*h2*h3*h5*(h0+h5)*(h2+h3)
         +(Pzp1*h5+Pzm1*h0)*h1*h2*h3*h4*(h1+h4)*(h2+h3)	
         )
      /((h0+h5)*(h1+h4)*(h2+h3)*(h0*h1*h4*h5+h0*h2*h3*h5+h1*h2*h3*h4));

   P[idx]=Csor*(tmp-P[idx])+P[idx];
   E1[idx]=(Pxp1-Pxm1)/(h2+h3);
   E2[idx]=(Pyp1-Pym1)/(h1+h4);
   E3[idx]=(Pzp1-Pzm1)/(h0+h5);
}

int Cylindrical::FindIdx(double tarx,double tary ,double tarz,int begin,int end)
{
   if(begin>=end)return XY::FindIdx(tarx,tary,begin,begin+x*y-1);
   int mid=((begin/(x*y)+end/(x*y))/2)*x*y;
   if(C3[mid]>=tarz)return FindIdx(tarx,tary,tarz,begin,mid);
   else return FindIdx(tarx,tary,tarz,mid+1,end);
}

double Cylindrical::GetData(double tarx, double tary, double tarz,int thing)
{
   int idx=FindIdx(tarx,tary,tarz,0,n);
   double ab=(tarx-C1[idx])/StepNext[idx];
   double aa=1-ab;
   double ba=(tary-C2[idx])/StepRight[idx];
   double bb=1-ba;
   double ac=(tarz-C3[idx])/StepUp[idx];
   double ca=1-ac;
   double tar0,tar1,tar2,tar3,tar4,tar5,tar6,tar7,*tar=NULL;
   switch(thing)
   {
      case 0:tar= Impurity;break;
      case 1:tar= P;break;
      case 2:tar= E1;break;
      case 3:tar= E2;break;
      case 4:tar=E3;break;
   }
   tar3=-1;
   tar5=-1;
   tar6=-1;
   tar7=-1;
   tar0=tar[idx];
   if(idx>=(n-x*y)){tar4=0;tar5=0;tar6=0;tar7=0;}
   else{tar4=tar[idx+x*y];}
   if(idx%(x*y)%x==x-1){tar2=0;tar3=0;tar6=0;tar7=0;}
   else{tar2=tar[idx+x];}
   if(idx%(x*y)/x==y-1){tar1=0;tar3=0;tar5=0;tar7=0;}
   if(tar3==-1)tar3=tar[idx+x+1];
   if(tar5==-1)tar5=tar[idx+x*y+1];
   if(tar6==-1)tar6=tar[idx+x*y+x];
   if(tar7==-1)tar7=tar[idx+x*y+x+1];
   return ((tar0*aa+tar1*ab)*ba+(tar2*aa+tar3*ab)*bb)*ac+((tar4*aa+tar5*ab)*ba+(tar6*aa+tar7*ab)*bb)*ca;
}

void Cylindrical::Save(const char * fout)
{
   XY::Save(fout);
   TFile *file=new TFile(fout,"update");
   TVectorD  v=*(TVectorD*)file->Get("v");
   v[9]=(double)z;
   v.Write();
   TTree * tree=(TTree*)file->Get("t");
   double E3s,C3s;
   tree->Branch("e3",&E3s,"e3/D"); // Electric field in z
   tree->Branch("c3",&C3s,"c3/D"); // persition in z
   for(int i=0;i<n;i++) {
      E3s=E3[i];
      C3s=C3[i];
      tree->Fill();
   }
   file->Write();
   file->Close();
   delete file;

}

void Cylindrical::Load(const char * fin)
{
   XY::Load(fin);
   TFile *file=new TFile(fin);
   TVectorD *v1=(TVectorD*)file->Get("v");
   double * v=v1->GetMatrixArray();
   z		=(int)	v[9];

   TChain *t =new TChain("t");
   t->Add(fin);
   double fEz,fPz;
   t->SetBranchAddress("c3",&fPz);
   t->SetBranchAddress("e3",&fEz);

   E3=new double[n];
   C3=new double[n];

   for (int i=0;i<n;i++) {
      t->GetEntry(i);
      E3[i]=fEz;
      C3[i]=fPz;
   }
   file->Close();
   delete file;
}

void Cylindrical::SetImpurity(TF3 * Im)
{
   for(int i=n;i-->0;) {
      Impurity[i]=Im->Eval((double)C1[i],(double)C2[i],(double)C3[i]);
   }
}

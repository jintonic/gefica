#include "X.h"

using namespace std;
using namespace GEFICA;

void X::Create(double steplength)
{
   E0=8.854;ER=16;Csor=1;Xlimit=0.0000001;

   E1=new double[n];
   C1=new double[n];
   P=new double[n];
   isbegin=new bool[n];
   StepNext=new double[n];
   StepBefore=new double[n];
   Impurity=new double[n];
   for (int i=n;i-->0;)
   {
      isbegin[i]=false;
      E1[i]=0;
      C1[i]=i*steplength;
      P[i]=0;
      StepNext[i]=steplength;
      StepBefore[i]=steplength;
      Impurity[i]=0;
   }
}

bool X::Iterate()
{
   int cnt=0,looptime=MaxIterations;
   while (cnt++<looptime)
   {
      XUpSum=0;
      XDownSum=0;
      for (int i=0;i<n;i++)
      {
         double tmp=P[i];
         Update(i);
         if(tmp>0)XDownSum+=tmp;
         else XDownSum=XDownSum-tmp;
         if(P[i]-tmp>0)XUpSum=XUpSum+P[i]-tmp;
         else XUpSum=XUpSum+tmp-P[i];
      }
      if(cnt%1000==0)cout<<cnt<<"  "<<XUpSum/XDownSum<<endl;
      if (XUpSum/XDownSum<Xlimit){return true;}
   }
   cout<<cnt<<endl;
   return false;
}

void X::Update(int idx)
{
   if (isbegin[idx])return;
   double density=Impurity[idx]*1.6/1000000000000;
   double h2=StepBefore[idx];
   double h3=StepNext[idx];
   double tmp=density/(E0*ER)*h2*h3+(h3*P[idx-1]+h2*P[idx+1])/(h2+h3);
   P[idx]=Csor*(tmp-P[idx])+P[idx];
   E1[idx]=(P[idx+1]-P[idx-1])/(h2+h3);
}

int X::FindIdx(double tarx,int begin,int end)
{
   if (begin>=end)return begin;
   int mid=(begin+end)/2;
   if(C1[mid]>=tarx)return FindIdx(tarx,begin,mid);
   else return FindIdx(tarx,mid+1,end);
}

double X::GetData(double tarx,int thing)
{
   int idx=FindIdx(tarx,0,n-1);
   if (idx==n)
   {
      switch (thing)
      {
         case 0:return Impurity[idx];
         case 1:return E1[idx];
         case 2:return P[idx];
      }
   }
   double ab=(tarx-C1[idx])/StepNext[idx];
   double aa=1-ab;
   switch(thing)
   {
      case 2:return E1[idx]*ab+E1[idx+1]*aa;
      case 1:return P[idx]*ab+C1[idx+1]*aa;
      case 0:return Impurity[idx]*ab+Impurity[idx+1]*aa;
   }
   return -1;
}
void X::Save(const char * fout)
{
   TFile * file=new TFile(fout,"recreate","data");
   TTree * tree=new TTree("t","1D");

   TVectorD v(10);

   v[7]=(double)x;
   v[8]=1;;
   v[9]=1;
   v[0]=(double)MaxIterations;
   v[1]=(double)n;
   v[2]=Csor;
   v[3]=E0;
   v[4]=ER;
   v[5]=XUpSum;
   v[6]=XDownSum;
   v.Write("v");
   bool isbegins;

   double E1s,C1s,Ps,StepNexts,StepBefores,impuritys;
   tree->Branch("e1",&E1s,"e1/D"); // Electric X in x
   tree->Branch("c1",&C1s,"c1/D"); // persition in x
   tree->Branch("p",&Ps,"p/D"); // electric potential
   tree->Branch("sn",&StepNexts,"StepNext/D"); // Step length to next point in x
   tree->Branch("sb",&StepBefores,"StepBefore/D"); // Step length to before point in x
   tree->Branch("ib",&isbegins,"isbegin/O"); // check is initial point
   tree->Branch("im",&impuritys,"impurity/D"); // Impurity
   for(int i=0;i<n;i++)
   {
      impuritys=Impurity[i];
      E1s=E1[i];
      C1s=C1[i];
      Ps=P[i];
      StepNexts=StepNext[i];
      StepBefores=StepBefore[i];
      tree->Fill();
   }
   file->Write();
   file->Close();
   delete file;
}
void X::Load(const char * fin)
{
   TFile *file=new TFile(fin);
   TVectorD *v1=(TVectorD*)file->Get("v");
   double * v=v1->GetMatrixArray();
   x		=(int)	v[7];
   MaxIterations	=(int)	v[0];
   n		=(int)	v[1];
   Csor		=	v[2];
   E0		=	v[3];
   ER		=	v[4];
   XUpSum	=	v[5];
   XDownSum	=	v[6];


   TChain *t =new TChain("t");
   t->Add(fin);
   bool fisbegin;
   double fE1,fC1,fP,fStepNext,fStepBefore,fimpurity;
   t->SetBranchAddress("c1",&fC1);
   t->SetBranchAddress("p",&fP);
   t->SetBranchAddress("sn",&fStepNext);
   t->SetBranchAddress("sb",&fStepBefore);
   t->SetBranchAddress("e1",&fE1);
   t->SetBranchAddress("ib",&fisbegin);
   t->SetBranchAddress("im",&fimpurity);

   E1=new double[n];
   C1=new double[n];
   P=new double[n];
   isbegin=new bool[n];
   StepNext=new double[n];
   StepBefore=new double[n];
   Impurity=new double[n];

   for (int i=0;i<n;i++)
   {
      t->GetEntry(i);
      E1[i]=fE1;
      C1[i]=fC1;
      P[i]=fP;
      isbegin[i]=fisbegin;  
      StepNext[i]=fStepNext;
      StepBefore[i]=fStepBefore;
      Impurity[i]=fimpurity;
   }
   file->Close();
   delete file;
}

void X::SetImpurity(TF1 * Im)
{
   for(int i=n;i-->0;)
   {
      Impurity[i]=Im->Eval((double)C1[i]);
   }
}

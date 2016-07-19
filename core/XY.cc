#include "XY.h"

using namespace GEFICA;

void XY::Create(double steplength)
{
  Field::Create(steplength);
  E2=new double[n];
  C2=new double[n];
  StepLeft=new double[n];
  StepRight=new double[n];
  for (int i=0;i<n;i++)
  {
    if(i>x-1)C2[i]=C2[i-x]+steplength;
    else C2[i]=0;
    if(i%x==0)C1[i]=0;
    else C1[i]=C1[i-1]+steplength;

    E2[i]=0;
    StepLeft[i]=steplength;
    StepRight[i]=steplength;
  }
}


void XY::Update(int idx)
{//need update
  if (isbegin[idx])return;
  double density=Impurity[idx]*1.6e12;
  double h2=StepBefore[idx];
  double h3=StepNext[idx];
  double h4=StepLeft[idx];
  double h1=StepRight[idx];
  double Pym1,Pyp1,Pxm1,Pxp1;
  if(idx>=x)Pym1=P[idx-x];
  else Pym1=P[idx];
  if(idx>=n-x)Pyp1=P[idx];
  else Pyp1=P[idx+x];
  if(idx%x==0)Pxm1=P[idx];
  else Pxm1=P[idx-1];
  if(idx%x==x-1)Pxp1=P[idx];
  else Pxp1=P[idx+1];
  double tmp=((h1+h4)*(h1*h2*h4*Pxp1+h1*h3*h4*Pxm1)+(h2+h3)*(h1*h2*h3*Pyp1+h2*h3*h4*Pym1)-density/(E0*ER)*(h1+h4)*(h2+h3)*h1*h2*h3*h4)/((h1+h4)*(h1*h2*h4+h1*h3*h4)+(h2+h3)*(h1*h2*h3+h2*h3*h4));
  P[idx]=Csor*(tmp-P[idx])+P[idx];
  E1[idx]=(Pxp1-Pxm1)/(h2+h3);
  E2[idx]=(Pyp1-Pym1)/(h1+h4);
}

int XY::FindIdx(double tarx,double tary ,int ybegin,int yend)
{
  if(ybegin>=yend)return X::FindIdx(tarx,ybegin,ybegin+x-1);
  int mid=((ybegin/x+yend/x)/2)*x;
  if(C2[mid]>=tary)return FindIdx(tarx,tary,ybegin,mid);
  else return FindIdx(tarx,tary,mid+1,yend);
}

double XY::GetData(double tarx, double tary, int thing)
{
  int idx=FindIdx(tarx,tary,0,n);
  double ab=(tarx-C1[idx])/StepNext[idx];
  double aa=1-ab;
  double ba=(tary-C2[idx])/StepRight[idx];
  double bb=1-ba;
  double tar0,tar1,tar2,tar3,*tar=NULL;
  switch(thing)
  {
    case 0:tar= Impurity;break;
    case 1:tar= P;break;
    case 2:tar= E1;break;
    case 3:tar= E2;break;
  }
  tar3=-1;
  tar0=tar[idx];
  if(idx/x+1==x){tar1=0;tar3=0;}
  else {tar1=tar[idx+1];}
  if(idx>n-x){tar2=0;tar3=0;}
  else {tar2=tar[idx+x];}
  if (tar3==-1)tar3=tar[idx+x+1];
  return (tar0*aa+tar1*ab)*ba+(tar2*aa+tar3*ab)*bb;
}
void XY::Save(const char * fout)
{
  Field::Save(fout);
  TFile *file=new TFile(fout,"update");
  TVectorD  v=*(TVectorD*)file->Get("v");
  v[8]=(double)y;
  v.Write();
  TTree * tree=(TTree*)file->Get("t");
  double E2s,C2s;
  tree->Branch("e2",&E2s,"e2/D"); // Electric field in y
  tree->Branch("c2",&C2s,"c2/D"); // persition in y
  for(int i=0;i<n;i++)
  {
    E2s=E2[i];
    C2s=C2[i];
    tree->Fill();
  }
  file->Write();
  file->Close();
  delete file;

}
void XY::Load(const char * fin)
{
  Field::Load(fin);
  TFile *file=new TFile(fin);
  TVectorD *v1=(TVectorD*)file->Get("v");
  double * v=v1->GetMatrixArray();
  y		=(int)	v[8];

  TChain *t =new TChain("t");
  t->Add(fin);
  double fEy,fPy;
  t->SetBranchAddress("c2",&fPy);
  t->SetBranchAddress("e2",&fEy);

  E2=new double[n];
  C2=new double[n];

  for (int i=0;i<n;i++)
  {
    t->GetEntry(i);
    E2[i]=fEy;
    C2[i]=fPy;
  }
  file->Close();
  delete file;
}
void XY::SetImpurity(TF2 * Im)
{
  for(int i=n;i-->0;)
  {
    Impurity[i]=Im->Eval((double)C1[i],(double)C2[i]);
  }
}

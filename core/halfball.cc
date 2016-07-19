#include "halfball.h"

using namespace GEFICA;

void halfball::Create(double steplength)
{
  n=x*y*z;
  Field2D::Create(steplength);
  E3=new double[n];
  C3=new double[n];
  StepUp=new double[n];
  StepDown=new double[n];
  for (int i=0;i<n;i++)
  {
    if(i/(x*y)==0)C3[i]=-3.1415926;
    else C3[i]=C3[i-x*y]+3.14159265*2/z;
    if((i%(x*y))/x!=0)C2[i]=C2[i-x]-3.1415926/(2*y-1);
    else C2[i]=1.570796;
    if(i%x==0)C1[i]=0;
    else C1[i]=C1[i-1]+steplength;

    E3[i]=0;
    StepLeft[i]=3.1415/z;
    StepRight[i]=3.1415/z;
    StepUp[i]=3.1415/(y*2-1);
    StepDown[i]=3.1415/(y*2-1);
  }
}


void halfball::Update(int idx)
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
  if(idx<x*y)Pzm1=P[idx+n-x*y];
  else Pzm1=P[idx-x*y];
  if(idx>=n-x*y)Pzp1=P[idx-(n-x*y)];
  else Pzp1=P[idx+x*y];
  if(idx%(x*y)>(x*y)-x-1)
  {
    if(idx<n/2)Pyp1=P[idx+n/2];
    else Pyp1=P[idx-n/2];
  }
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

int halfball::FindIdx(double tarx,double tary ,double tarz,int begin,int end)
{
  if(begin>=end)return Field2D::FindIdx(tarx,tary,begin,begin+x*y-1);
  int mid=((begin/(x*y)+end/(x*y))/2)*x*y;
  if(C3[mid]>=tarz)return FindIdx(tarx,tary,tarz,begin,mid);
  else return FindIdx(tarx,tary,tarz,mid+1,end);
}

double halfball::GetData(double tarx, double tary, double tarz,int thing)
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
  if(tary==0)return (tar[x*y-1]+tar[x*y-1+n/2])/2;
  tar3=-1;
  tar5=-1;
  tar6=-1;
  tar7=-1;
  tar0=tar[idx];
  if(idx>=(n-x*y)){tar4=tar[idx-n+x*y];tar5=tar[idx-n+x*y+1];tar6=tar[idx-n+x*y+x];tar7=tar[idx-n+x*y+x+1];}
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
void halfball::Save(const char * fout)
{
  Field2D::Save(fout);
  TFile *file=new TFile(fout,"update");
  TVectorD  v=*(TVectorD*)file->Get("v");
  v[9]=(double)z;
  v.Write();
  TTree * tree=(TTree*)file->Get("t");
  double E3s,C3s;
  tree->Branch("e3",&E3s,"e3/D"); // Electric field in z
  tree->Branch("c3",&C3s,"c3/D"); // persition in z
  for(int i=0;i<n;i++)
  {
    E3s=E3[i];
    C3s=C3[i];
    tree->Fill();
  }
  file->Write();
  file->Close();
  delete file;

}
void halfball::Load(const char * fin)
{
  Field2D::Load(fin);
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

  for (int i=0;i<n;i++)
  {
    t->GetEntry(i);
    E3[i]=fEz;
    C3[i]=fPz;
  }
  file->Close();
  delete file;
}
void halfball::SetImpurity(TF3 * Im)
{
  for(int i=n;i-->0;)
  {
    Impurity[i]=Im->Eval((double)C1[i],(double)C2[i],(double)C3[i]);
  }
}

#include <Field.h>

void Field::Create(double steplength)
{
  Ex=new double[n];
  Px=new double[n];
  f=new double[n];
  isbegin=new bool[n];
  StepNext=new double[n];
  StepBefore=new double[n];
  for (int i=n;i-->0;)
  {
    isbegin[i]=false;
    Ex[i]=0;
    Px[i]=i*steplength;
    f[i]=0;
    StepNext[i]=steplength;
    StepBefore[i]=steplength;
  }
}

bool Field::Iteration()
{
  while (cnt<looptime)
  {
    XUpSum=0;
    XDownSum=0;
    for (int i=0;i<n;i++)
    {
      double tmp=f[i];
      update(i);
      if(tmp>0)XDownSum+=tmp;
      else XDownSum=XDownSum-tmp;
      if(f[i]-tmp>0)XUpSum=XUpSum+f[i]-tmp;
      else XUpSum=XUpSum+tmp-f[i];
    } 
    if (XUpSum/XDownSum<Xlimit)return true;
  }return false;
}

void Field::Update(int idx)
{
  if (isbegin[idx])return;
  double density=16;
  double h2=StepBefore[idx];
  double h3=StepNext[idx];
  double tmp=density/(E0*ER)*h2*h3+(h3*f[idx-1]+h2*f[idx+1])/(h2+h3);
  f[idx]=Csor*(tmp-f[idx])+f[idx];
  Ex[idx]=(f[idx+1]-f[idx-1])/(h2+h3);
}

int Field::FindIdx(double tarx,int begin,int end)
{
  if (begin>=end)return begin;
  int mid=(begin+end)/2;
  if(Px[mid]>=tarx)return FindIdx(tarx,begin,mid);
  else return FindIdx(tarx,mid+1,end);
}

double Field::GetData(double x,int thing)
{
  int idx=FindIdx(x,0,x);
  double ab=(x-Px[idx])/StepNext[idx];
  double aa=1-ab;
  switch(thing)
  {
    case 1:return Ex[idx]*ab+Ex[idx+1]*aa;
    case 2:return Px[idx]*ab+Px[idx+1]*aa;
    case 3:return StepNext[idx+1]-x;
    case 4:return x-StepBefore[idx];
  }
  return -1;
}
void Field::save(const char * fout)
{
  TFile * file=new TFile(fout,"recreate","data");
  TTree * tree=new TTree("t","1D");
  TVecotorD *v=new TVectorD(8);
  v[0]=(double)x;
  v[1]=(double)tlimit;
  v[2]=(double)n;
  v[3]=Csor;
  v[4]=E0;
  v[5]=ER;
  v[6]=XUpSum;
  v[7]=XDownSum;
  v->write("v");
  bool isbegins;
  double Exs,Ps,fs,StepNexts,StepBefores;
  tree->Branch("Ex",&Exs,"Exs/D");
  tree->Branch("Px",&Ps,"Ps/D"); 
  tree->Branch("f",&fs,"fs/D");
  tree->Branch("StepNext",&StepNexts,"StepNexts/D");
  tree->Branch("StepBefore",&StepBefores,"StepBefores/D");
  tree->Branch("isbegin",&isbegins,"isbegins/O");
  for(int i=0;i<x;i++)
  {
    Exs=Ex[i];
    Ps=Px[i];
    fs=f[i];
    StepNexts=StepNext[i];
    StepBefores=StepBefore[i];
    tree.fill();
  }
  file.write;
  file.close();
  delete file;
}
void Field:load(const char * fin)
{
  TFile *file=new TFile(fin);
  TVector *v=(TVectorD*)file->get("v");
  x		=(int)	v[0];
  tlimit	=(int)	v[1];
  n		=(int)	v[2];
  Csor		=	v[3];
  E0		=	v[4];
  ER		=	v[5];
  XUpSum	=	v[6];
  XDownSum	=	v[7];

  TChain *t =new TChain("t");
  t->Add(fin);
  bool fisbegin;
  double fEx,fPx,ff,fStepNext,fStepBefore;
  t->SetBranchAddress("Px",&fP);
  t->SetBranchAddress("f",&ff);
  t->SetBranchAddress("StepNext",&fStepNext);
  t->SetBranchAddress("StepBefore",&fStepBefore);
  t->SetBranchAddress("Ex",&fEx);
  t->SetBranchAddress("isbegin",&fisbegin);
  t->GetEntry(0);
  Ex=new double[n];
  Px=new double[n];
  f=new double[n];
  isbegin=new bool[n];
  StepNext=new double[n];
  StepBefore=new double[n];

  for (int i=0;i<n;i++)
  {
    Ex[i]=fEx;
    Px[i]=fPx;
    f[i]=ff;
    isbegin[i]=fisbegin;  
    StepNext[i]=fStepNext;
    StepBefore[i]=fStepBefore;
  }
}

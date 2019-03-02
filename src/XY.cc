#include <TF3.h>
#include <TFile.h>
#include <TChain.h>
#include <TVectorD.h>

#include "XY.h"
#include "Units.h"
using namespace GeFiCa;

XY::XY(int nx, int ny, const char *name, const char *title)
   : X(nx*ny, name, title), n2(ny)
{
   n1=nx; // n1 is set to nx*ny through X constructor, it is fixed here
}
//_____________________________________________________________________________
//
XY::~XY()
{
   if (fE2) delete[] fE2;
   if (fC2) delete[] fC2;
   if (fdC2m) delete[] fdC2m;
   if (fdC2p) delete[] fdC2p;
}
//_____________________________________________________________________________
//
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
      fdC2m[i]=steplength2;
      fdC2p[i]=steplength2;
   }
}
//_____________________________________________________________________________
//
#include <iostream>
using namespace std;
void XY::DoSOR2(int idx)
{
   // 2nd-order Runge-Kutta Successive Over-Relaxation
   if (fIsFixed[idx])return;
   double density=fImpurity[idx]*Qe;
   double h2=fdC1m[idx];
   double h3=fdC1p[idx];
   double h4=fdC2m[idx];
   double h1=fdC2p[idx];
   double pym,pyp,pxm,pxp;
   if(idx>=n1)pym=fV[idx-n1];
   else pym=fV[idx];
   if(idx>=n-n1)pyp=fV[idx];
   else pyp=fV[idx+n1];
   if(idx%n1==0)pxm=fV[idx];
   else pxm=fV[idx-1];
   if(idx%n1==n1-1)pxp=fV[idx];
   else pxp=fV[idx+1];
   double tmp=(density/epsilon+(pxp/h2+pxm/h3)*2/(h2+h3)+(pyp/h1+pym/h4)*2/(h1+h4))/
      ((1/h2+1/h3)*2/(h2+h3)+(1/h1+1/h4)*2/(h1+h4));
   //find minmium and maxnium of all five grid, the new one should not go overthem.
   //find min
   double min=pxm;
   double max=pxm;
   if(min>pxp)min=pxp;
   if (min>pyp)min=pyp;
   if (min>pym)min=pym;
   
   //find max
   if(max<pxp)max=pxp;
   if (max<pyp)max=pyp;
   if (max<pym)max=pym;
   //if tmp is greater or smaller than max and min, set tmp to it.
   
   //fV[idx]=Csor*(tmp-fV[idx])+fV[idx];
   //if need calculate depleted voltage
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
int XY::FindIdx(double tarx,double tary ,int ybegin,int yend)
{
   //search using binary search
   // if(ybegin>=yend)cout<<"to x"<<ybegin<<" "<<yend<<endl;;
   if(ybegin>=yend)return X::FindIdx(tarx,yend*n1,(yend+1)*n1-1);
   int mid=((ybegin+yend)/2);
   if(fC2[mid*n1]>=tary){//cout<<"firsthalf"<<ybegin<<" "<<yend<<endl; 
      return FindIdx(tarx,tary,ybegin,mid);
   }
   else{//cout<<"senondhalf"<<ybegin<<" "<<yend<<endl; 
      return FindIdx(tarx,tary,mid+1,yend);}
}
//_____________________________________________________________________________
//
double XY::GetData(double tarx, double tary, EOutput output)
{
   int idx=FindIdx(tarx,tary,0,n2-1);

   //test
   //cout<<"index:"<<idx<<endl;
   //cout<<"(0,0)c1: "<<fC1[idx]<<" c2: "<<fC2[idx]<<" p: "<<fV[idx]<<endl;
   //cout<<"(0,1)c1: "<<fC1[idx-1]<<" c2: "<<fC2[idx-1]<<" p: "<<fV[idx-1]<<endl;
   //cout<<"(1,0)c1: "<<fC1[idx-n1]<<" c2: "<<fC2[idx-n1]<<" p: "<<fV[idx-n1]<<endl;
   //cout<<"(1,1)c1: "<<fC1[idx-n1-1]<<" c2: "<<fC2[idx-n1-1]<<" p: "<<fV[idx-n1-1]<<endl;
   //
   //cout<<idx<<" "<<n<<endl;
   double ab=(-tarx+fC1[idx])/fdC1m[idx];
   double aa=1-ab;
   double ba=(-tary+fC2[idx])/fdC2m[idx];
   //cout<<"left"<<fdC2m[idx]<<endl;
   //cout<<"right "<<fdC2p[idx]<<endl;
   //cout<<"next"<<fdC1p[idx]<<endl;
   double bb=1-ba;
   double tar0,tar1,tar2,tar3,*tar=NULL;
   switch(output) {
      case 0:tar= fImpurity;break;
      case 1:tar= fV;break;
      case 2:tar= fE1;break;
      case 3:tar= fE2;break;
      default:break;
   }
   tar3=-1;
   tar0=tar[idx];
   if (idx%n1==0){tar1=0;tar3=0;}
   else {tar1=tar[idx-1];}
   if(idx<n1) {tar2=0;tar3=0;}
   else {tar2=tar[idx-n1];}
   if (tar3==-1)tar3=tar[idx-n1-1];
   //cout<<tar0<<" "<<tar1<<" "<<tar2<<" "<<tar3<<endl;
   //cout<<tarx<<", "<<tary<<endl;
   //cout<<aa<<" "<<ab<<" "<<ba<<" "<<bb<<endl;
   //
   //if (fC1[idx]>0. && fC2[idx]>3.445){
   //   cout<<tary-0.21<<endl;
   //   abort();
   //}
   return (tar1*ab+tar0*aa)*bb+(tar3*ab+tar2*aa)*ba;
}
//_____________________________________________________________________________
//
bool XY::CalculateField(int idx)
{
   if (!X::CalculateField(idx)) return false;
   if (fdC2p[idx]==0 || fdC2m[idx]==0) return false;

   if (idx%(n1*n2)<n1) // C2 lower boundary
      fE2[idx]=(fV[idx]-fV[idx+n1])/fdC2p[idx];
   else if (idx%(n1*n2)>=n-n1) // C2 upper boundary
      fE2[idx]=(fV[idx]-fV[idx-n1])/fdC2m[idx];
   else { // bulk
      fE2[idx]=(fV[idx-n1]-fV[idx+n1])/(fdC2m[idx]+fdC2p[idx]);
   }
   return true;
}
//_____________________________________________________________________________
//
TTree* XY::GetTree(bool createNew)
{
   X::GetTree(createNew); // create tree

   double e2,c2,dc2p,dc2m;
   TBranch *be2 = fTree->Branch("e2",&e2,"e2/D");
   TBranch *bc2 = fTree->Branch("c2",&c2,"c2/D");
   TBranch *bp2 = fTree->Branch("dc2p",&dc2p,"dc2p/D");
   TBranch *bm2 = fTree->Branch("dc2m",&dc2m,"dc2m/D");

   for(int j=0;j<n;j++) {
      e2= fE2[j];
      c2= fC2[j];
      dc2p=fdC2p[j];
      dc2m=fdC2m[j];

      be2->Fill();
      bc2->Fill();
      bp2->Fill();
      bm2->Fill();
   }

   fTree->ResetBranchAddresses(); // disconnect from local variables
   return fTree;
}

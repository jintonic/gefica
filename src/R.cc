#include <iostream>
using namespace std;

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TVectorD.h>
#include <TF1.h>

#include "R.h"
#include "Units.h"
using namespace GeFiCa;

void R::DoSOR2(int idx)
{
   if (fIsFixed[idx])return ;
   double density=fImpurity[idx]*Qe;
   double h2=fdC1m[idx];
   double h3=fdC1p[idx];
   double p2=fV[idx-1];
   double p3=fV[idx+1];
   double tmp=(+density/epsilon*(h2+h3)*0.5+1/fC1[idx]*(p3-p2)
         +p3/h2+p2/h3)/(1/h2+1/h3);

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
   } else fIsDepleted[idx]=true;

   if(fIsDepleted[idx]||V0==V1) fV[idx]=tmp;
}
//_____________________________________________________________________________
//
void R::DoSOR4(int idx)
{ 
   if (fIsFixed[idx])return;

   double density=fImpurity[idx]*1.6e-19;
   double h2=fdC1m[idx];
   double h3=fdC1p[idx];
   double h1=h2;
   double xp2,xm2,xm1,xp1;
   xm1=fV[idx-1];
   xp1=fV[idx+1];
   if(idx>1)xm2=fV[idx-2];
   else {DoSOR2(idx);return; } 
   if(idx<n-2)xp2=fV[idx+2];
   else {DoSOR2(idx);return;}
   double tmp=(-1/12*xp2+4/3*xp1+4/3*xm1-1/12*xm2-density/epsilon*h1*h1)*2/5;
   fV[idx]=Csor*(tmp-fV[idx])+fV[idx];

   fE1[idx]=(fV[idx+1]-fV[idx-1])/(h2+h3);
}


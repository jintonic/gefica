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

void R::SOR2(int idx,bool elec)
{
   if (fIsFixed[idx])return ;
   double density=fImpurity[idx]*Qe;
   double h2=fdC1m[idx];
   double h3=fdC1p[idx];
   double tmp=(+density/epsilon*(h2+h3)*0.5+1/fC1[idx]*(fPotential[idx+1]-fPotential[idx-1])
         +fPotential[idx+1]/h2+fPotential[idx-1]/h3)/(1/h2+1/h3);

   //fPotential[idx]=Csor*(tmp-fPotential[idx])+fPotential[idx];
   double oldP=fPotential[idx];
      tmp=Csor*(tmp-oldP)+oldP;
   
   if(tmp<min)
   {
      fPotential[idx]=min;
      fIsDepleted[idx]=false;
   }
   else if(tmp>max)
   {
      fPotential[idx]=max;
      fIsDepleted[idx]=false;
   }
   else
      fIsDepleted[idx]=true;
   if(fIsDepleted[idx]||!NotImpurityPotential)
   {
      //over relax
      fPotential[idx]=tmp;
   }
}
//_____________________________________________________________________________
//
void R::SOR4(int idx)
{ 
  if (fIsFixed[idx])return;

   double density=fImpurity[idx]*1.6e-19;
   double h2=fdC1m[idx];
   double h3=fdC1p[idx];
   double h1=h2;
   double xp2,xm2,xm1,xp1;
   xm1=fPotential[idx-1];
   xp1=fPotential[idx+1];
   if(idx>1)xm2=fPotential[idx-2];
   else {SOR2(idx,0);return; } 
   if(idx<n-2)xp2=fPotential[idx+2];
   else {SOR2(idx,0);return;}
   double tmp=(-1/12*xp2+4/3*xp1+4/3*xm1-1/12*xm2-density/epsilon*h1*h1)*2/5;
   fPotential[idx]=Csor*(tmp-fPotential[idx])+fPotential[idx];

   fE1[idx]=(fPotential[idx+1]-fPotential[idx-1])/(h2+h3);
}


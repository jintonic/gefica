#include <iostream>
using namespace std;

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TVectorD.h>
#include <TF1.h>

#include "Rho.h"
#include "Units.h"
using namespace GeFiCa;

void Rho::SOR2(int idx,bool elec)
{
   if (fIsFixed[idx])return ;
   double density=fImpurity[idx]*Qe;
   double h2=fdC1m[idx];
   double h3=fdC1p[idx];
   double p2=fPotential[idx-1];
   double p3=fPotential[idx+1];
   //double tmp=-density/epsilon*h2*h3/2+(fPotential[idx-1]-fPotential[idx+1])/fC1[idx]*h2*h3/(h2+h3)+(h3*fPotential[idx-1]+h2*fPotential[idx+1])/(h2+h3);
   double tmp=(+density/epsilon*(h2+h3)*0.5+0.5/fC1[idx]*(p3-p2)
         +p3/h2+p2/h3)/(1/h2+1/h3);
    // over-relaxation if Csor>1
   //find minmium and maxnium of all five grid, the new one should not go overthem.
   //find min
   double min=p2;
   double max=p2;
   if(min>p3)min=p3;
   //find max
   if(max<p3)max=p3;
//if tmp is greater or smaller than max and min, set tmp to it.
   
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
void Rho::SOR4(int idx)
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


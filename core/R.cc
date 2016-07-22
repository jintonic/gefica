#include <iostream>
using namespace std;

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TVectorD.h>
#include <TF1.h>

#include "R.h"
using namespace GEFICA;


void R::RK2(int idx)
{
   if (fIsFixed[idx])return ;
   double density=fImpurity[idx]*1.6e-19;
   double h2=fDistanceToPrevious[idx];
   double h3=fDistanceToNext[idx];
   double tmp=-density/epsilon*h2*h3/2-(fPotential[idx-1]-fPotential[idx+1])/fC1[idx]*h2*h3/(h2+h3)+(h3*fPotential[idx-1]+h2*fPotential[idx+1])/(h2+h3);
   // over-relaxation if Csor>1
   fPotential[idx]=Csor*(tmp-fPotential[idx])+fPotential[idx];

   fE1[idx]=(fPotential[idx+1]-fPotential[idx-1])/(h2+h3);
}

void R::RK4(int idx)
{ 
  if (fIsFixed[idx])return;

   double density=fImpurity[idx]*1.6e-19;
   double h2=fDistanceToPrevious[idx];
   double h3=fDistanceToNext[idx];
   double h1=h2;
   double h4=h3;
   double xp2,xm2,xm1,xp1;
   xm1=fPotential[idx-1];
   xp1=fPotential[idx+1];
   if(idx>1)xm2=fPotential[idx-2];
   else {RK2(idx);return; } 
   if(idx<n-2)xp2=fPotential[idx+2];
   else {RK2(idx);return;}
   double tmp=(-1/12*xp2+4/3*xp1+4/3*xm1-1/12*xm2-density/epsilon*h1*h1)*2/5;
   fPotential[idx]=Csor*(tmp-fPotential[idx])+fPotential[idx];

   fE1[idx]=(fPotential[idx+1]-fPotential[idx-1])/(h2+h3);
}


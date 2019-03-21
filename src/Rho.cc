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

void Rho::OverRelaxAt(int idx)
{
   if (fIsFixed[idx])return ;
   double density=fImpurity[idx]*Qe;
   double h2=dC1m[idx];
   double h3=dC1p[idx];
   double p2=V[idx-1];
   double p3=V[idx+1];
   double tmp=(+density/epsilon*(h2+h3)*0.5+0.5/C1[idx]*(p3-p2)
         +p3/h2+p2/h3)/(1/h2+1/h3);
   // over-relaxation if RelaxationFactor>1
   //find minmium and maxnium of all five grid, the new one should not go overthem.
   //find min
   double min=p2;
   double max=p2;
   if(min>p3)min=p3;
   //find max
   if(max<p3)max=p3;
   //if tmp is greater or smaller than max and min, set tmp to it.

   //V[idx]=RelaxationFactor*(tmp-V[idx])+V[idx];
   double oldP=V[idx];
   tmp=RelaxationFactor*(tmp-oldP)+oldP;

   if(tmp<min) {
      V[idx]=min;
      fIsDepleted[idx]=false;
   } else if(tmp>max) {
      V[idx]=max;
      fIsDepleted[idx]=false;
   } else
      fIsDepleted[idx]=true;

   if(fIsDepleted[idx]||Bias[0]==Bias[1]) V[idx]=tmp;
}
//_____________________________________________________________________________
//
void Rho::DoSOR4(int idx)
{ 
   if (fIsFixed[idx])return;

   double density=fImpurity[idx]*1.6e-19;
   double h2=dC1m[idx];
   double h3=dC1p[idx];
   double h1=h2;
   double xp2,xm2,xm1,xp1;
   xm1=V[idx-1];
   xp1=V[idx+1];
   if(idx>1)xm2=V[idx-2];
   else {OverRelaxAt(idx);return; } 
   if(idx<fN-2)xp2=V[idx+2];
   else {OverRelaxAt(idx);return;}
   double tmp=(-1/12*xp2+4/3*xp1+4/3*xm1-1/12*xm2-density/epsilon*h1*h1)*2/5;
   V[idx]=RelaxationFactor*(tmp-V[idx])+V[idx];

   E1[idx]=(V[idx+1]-V[idx-1])/(h2+h3);
}


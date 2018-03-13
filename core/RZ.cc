#include <TF2.h>
#include <TFile.h>
#include <TChain.h>
#include <TVectorD.h>

#include "RZ.h"
#include "Units.h"

using namespace GeFiCa;

//_____________________________________________________________________________
//
#include <iostream>
using namespace std;
void RZ::SOR2(int idx,bool elec)
{
   // 2nd-order Successive Over-Relaxation
   if (fIsFixed[idx])return;
   double density=fImpurity[idx]*Qe;
   double h2=fDistanceToPrevious[idx];
   double h3=fDistanceToNext[idx];
   double h4=fDistanceToLeft[idx];
   double h1=fDistanceToRight[idx];
   double Pym1,Pyp1,Pxm1,Pxp1;
   if(idx>=n1)Pym1=fPotential[idx-n1];
   else Pym1=fPotential[idx];
   if(idx>=n-n1)Pyp1=fPotential[idx];
   else Pyp1=fPotential[idx+n1];
   if(idx%n1==0)Pxm1=fPotential[idx];
   else Pxm1=fPotential[idx-1];
   if(idx%n1==n1-1)Pxp1=fPotential[idx];
   else Pxp1=fPotential[idx+1];
   double tmp=(density/epsilon+1/fC1[idx]*(Pxp1-Pxm1)/(h2+h3)+(Pxp1/h3+Pxm1/h2)*2/(h2+h3)+(Pyp1/h1+Pym1/h4)*2/(h1+h4))/
      ((1/h2+1/h3)*2/(h2+h3)+(1/h1+1/h4)*2/(h1+h4));
   fPotential[idx]=Csor*(tmp-fPotential[idx])+fPotential[idx];
   if(elec) {
      fE1[idx]=(Pxp1-Pxm1)/(h2+h3);
      fE2[idx]=(Pyp1-Pym1)/(h1+h4);
   }
}

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
   double density=fImpurity[idx]*Qe;
   double drm=fdC1m[idx]; // dr_minus
   double drp=fdC1p[idx];
   double dzm=fdC2m[idx];
   double dzp=fdC2p[idx];
   double pzm,pzp,prm,prp; // pzm: potential_z_plus
   if(idx>=n1)pzm=fPotential[idx-n1];
   else pzm=fPotential[idx];
   if(idx>=n-n1)pzp=fPotential[idx];
   else pzp=fPotential[idx+n1];
   if(idx%n1==0)prm=fPotential[idx];
   else prm=fPotential[idx-1];
   if(idx%n1==n1-1)prp=fPotential[idx];
   else prp=fPotential[idx+1];
   if (!fIsFixed[idx]) {
      double tmp=(density/epsilon
            + 1/fC1[idx]*(prp-prm)/(drm+drp) +(prp/drp+prm/drm)*2/(drm+drp)
            + (pzp/dzp+pzm/dzm)*2/(dzp+dzm))/
         ((1/drm+1/drp)*2/(drm+drp)+(1/dzp+1/dzm)*2/(dzp+dzm));
      fPotential[idx]=Csor*(tmp-fPotential[idx])+fPotential[idx];
   }
   if (elec) {
      fE1[idx]=(prm-prp)/(drm+drp);
      fE2[idx]=(pzm-pzp)/(dzp+dzm);

      if (idx<n1) fE2[idx]=(fPotential[idx]-pzp)/dzp; // lower border
      if (idx>=n-n1) fE2[idx]=(fPotential[idx]-pzm)/dzm; // upper border
      if (idx%n1==0) fE1[idx]=(fPotential[idx]-prp)/drp; // left border
      if (idx%n1==n1-1) fE1[idx]=(fPotential[idx]-prm)/drm; // right border
   }
}

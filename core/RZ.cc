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
   if (fIsFixed[idx])return; 
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
   double tmp=(density/epsilon
         + 1/fC1[idx]*(prp-prm)/(drm+drp) +(prp/drp+prm/drm)*2/(drm+drp)
         + (pzp/dzp+pzm/dzm)*2/(dzp+dzm))/
      ((1/drm+1/drp)*2/(drm+drp)+(1/dzp+1/dzm)*2/(dzp+dzm));
   //find minmium and maxnium of all five grid, the new one should not go overthem.
   //find min
   double min=prm;
   double max=prm;
   if(min>prp)min=prp;
   if (min>pzp)min=pzp;
   if (min>pzm)min=pzm;
   //find max
   if(max<prp)max=prp;
   if (max<pzp)max=pzp;
   if (max<pzm)max=pzm;
//if tmp is greater or smaller than max and min, set tmp to it.
   
      //over relax
   //fPotential[idx]=Csor*(tmp-fPotential[idx])+fPotential[idx];
   tmp=Csor*(tmp-fPotential[idx])+fPotential[idx];
   
   if(tmp<min)fPotential[idx]=min;
   else if(tmp>max)fPotential[idx]=max;
   else fPotential[idx]=tmp;



}


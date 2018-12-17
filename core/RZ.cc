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
void RZ::SOR2(int idx,bool NotImpurityPotential)
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
   //if need calculate depleted voltage
   double oldP=fPotential[idx];
   
   if(tmp<min)
   {
      fPotential[idx]=min;
      DepletedData[idx]=false;
   }
   else if(tmp>max)
   {
      fPotential[idx]=max;
      DepletedData[idx]=false;
   }
   else
      DepletedData[idx]=true;
   if(DepletedData[idx]||!NotImpurityPotential)
   {
      tmp=Csor*(tmp-oldP)+oldP;
      fPotential[idx]=tmp;
   }



}
void RZ::CalculateCapacitance()
{
   /*
   http://bgaowww.physics.utoledo.edu/teaching/LectureNotes/Phys2080/Chapter16.htm
   https://en.wikipedia.org/wiki/Electric_field#Energy_in_the_electric_field
   use Energy in a charged capacitor= total energy U stored in the electric field in a given volume V 
   to find capacitance
   */
   //only work when impurity are zero
   for(int i=0;i<n;i++)
   {
      if(fImpurity[i]!=0)
      {
         //impurity not clear,return
         return;
      }
   }
   double V=1*volt;
   //debug:cout<<V<<endl;
   double SumofElectricField=0;
   for(int i=0;i<n;i++)
   {
      if(fC1[i]<0)continue;
      //integral over electric field
      double e1=fE1[i];
      double e2=fE2[i];
      double dr=fdC1p[i];
      double dz=fdC2p[i];
      SumofElectricField+=(e1*e1+e2*e2)*fC1[i]*dr*dz;

   }
   double Capacitance=SumofElectricField*2*3.14159*epsilon/V/V;
   cout<<"Calculated capacitance is "<<Capacitance/pF<<endl;
}


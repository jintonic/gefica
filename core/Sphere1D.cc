#include <iostream>
using namespace std;

#include "Sphere1D.h"
using namespace GEFICA;
//____________________________________________________
//a Planar detector under 1D coordinate system with analyic
ClassImp(Sphere1D)
void Sphere1D::initialize()
{
   // The step length is calculated with the following equation:
   // BEGIN_HTML
   // <pre>
   //      double stepLength=(innerR-outterR)/(n-1);
   // </pre>
   // END_HTML
   // If the inner radius is not larger than the outer radius,
   // no grid will be created
   if (outterR>=innerR) {
      Warning("CreateGridWithFixedStepLength",
            "Lower bound (%f) >= upper bound (%f)! No grid is created!",
            innerR, outterR);
      return;
   }
   double steplength=(innerR-outterR)/(n-1);
   SetStepLength(steplength);
   for(int i=n;i-->0;)fC1[i]=fC1[i]+outterR;
   fIsFixed[0]=true;
   fIsFixed[n-1]=true;
   double slope = (cathode_voltage-annode_voltage)/(n-1);
   for (int i=0; i<n; i++) fPotential[i]=annode_voltage+slope*i;
}
#include  <cmath>
bool Sphere1D::Analyic()
{
   double density=fImpurity[1]*1.6e-19;
   double cnst1=(fPotential[n-1]-fPotential[0]+density/epsilon/6*(fC1[n-1]*fC1[n-1]-fC1[0]*fC1[0]))/(pow(fC1[n-1],-1)-pow(fC1[0],-1));
   cout<<fPotential[0]<<endl;
   double cnst2=fPotential[0]-density/epsilon/6*fC1[0]*fC1[0]-cnst1*pow(fC1[0],-1);//fPotential[0]-density*fC1[0]*fC1[0]/epsilon/4-cnst1*log(fC1[0]);
   for (int i=0; i<n; i++) {
      fPotential[i] = fImpurity[i]*1.6e-19/6/epsilon*fC1[i]*fC1[i]+cnst1*pow(fC1[i],-1)+cnst2;
      cout<<fPotential[i]<<endl;
      fE1[i]=(fPotential[i+1]-fPotential[i-1])/(fDistanceToNext[i]+fDistanceToPrevious[i]);
   }
   return true;
}
bool Sphere1D::CalculateField(EMethod method)
{
  if(!fIsLoaded)initialize();
  return R::CalculateField(method);
}

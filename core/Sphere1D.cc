#include <iostream>
using namespace std;

#include "Sphere1D.h"
using namespace GEFICA;

void Sphere1D::SetVoltage(double anode_voltage, double cathode_voltage)
{
   fIsFixed[0]=true;
   fIsFixed[n-1]=true;
   double slope = (cathode_voltage-anode_voltage)/(n-1);
   for (int i=0; i<n; i++) {
      fPotential[i]=anode_voltage+slope*i;
   }
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

void Sphere1D::CreateGridWithFixedStepLength()
{
   // The step length is calculated with the following equation:
   // BEGIN_HTML
   // <pre>
   //      double stepLength=(OuterRadius-InnerRadius)/(n-1);
   // </pre>
   // END_HTML
   // If the inner radius is not larger than the outer radius,
   // no grid will be created
   if (InnerRadius>=OuterRadius) {
      Warning("CreateGridWithFixedStepLength",
            "Inner radius (%f) >= outer radius (%f)! No grid is created!",
            InnerRadius, OuterRadius);
      return;
   }
   double stepLength=(OuterRadius-InnerRadius)/(n-1);
   X::CreateGridWithFixedStepLength(stepLength);
   for(int i=0;i<n;i++) fC1[i]=fC1[i]+r0;
}

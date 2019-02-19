#include "Units.h"
#include "Sphere3D.h"
using namespace GeFiCa;

Sphere3D::Sphere3D(int nr, int nt, int np, const char *name, const char *title)
   : RThetaPhi(nr, nt, np, name, title), OuterRadius(3*cm), InnerRadius(0.3*cm)
{};
//______________________________________________________________________________
//
void Sphere3D::Initialize()
{
   if (InnerRadius>=OuterRadius) Fatal("Initialize",
         "Inner R (%.1f) >= outer R (%.1f)! Abort!", InnerRadius, OuterRadius);

   double stepLength=(OuterRadius-InnerRadius)/(n1-1);
   SetStepLength(stepLength,3.14159265/n2,3.14159265*2/n3);
   for(int i=n;i-->0;) fC1[i]=fC1[i]+InnerRadius;
   for(int i=n;i-->0;) fC2[i]=fC2[i]+3.14159265/2/n2;

   for (int i=0; i<n; i+=n1) {
      fV[i]=V0;
      fV[i+n1-1]=V1;
      fIsFixed[i]=true;
      fIsFixed[i+n1-1]=true;
   }
}

#include "Units.h"
#include "SegmentedInZ.h"
using namespace GeFiCa;

SegmentedInZ::SegmentedInZ(int nr, int nt, const char *name, const char *title)
   : RhoZ(nr, nt, name, title), InnerR(0.5*cm), OuterR(3*cm), Z(3*cm), Nseg(3),
   SegmentId(1) {};
//______________________________________________________________________________
//
void SegmentedInZ::Initialize()
{
   if (InnerR>=OuterR) Fatal("Initialize",
            "Inner R (%.1f) >= outer R (%.1f)! Abort!", InnerR, OuterR);

   double stepLength=(OuterR-InnerR)/(n1-1);
   SetStepLength(stepLength,2*3.14159265/n2);

   for(int i=n;i-->0;) fC1[i]=fC1[i]+InnerR;
   for (int i=0;i<n;i=i+n1) {
      fV[i]=V1;
      fV[i+n1-1]=V0;
      fIsFixed[i]=true;
      fIsFixed[i+n1-1]=true;
   }
}

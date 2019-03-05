#include "Units.h"
#include "SegmentedInZ.h"
using namespace GeFiCa;

SegmentedInZ::SegmentedInZ(int nr, int nt, const char *name, const char *title)
   : RhoZ(nr, nt, name, title), InnerR(0.5*cm), OuterR(3*cm), Z(3*cm), Nseg(3),
   SegmentId(1) {};
//______________________________________________________________________________
//
void SegmentedInZ::InitializeGrid()
{
   if (InnerR>=OuterR) Fatal("InitializeGrid",
            "Inner R (%.1f) >= outer R (%.1f)! Abort!", InnerR, OuterR);

   double stepLength=(OuterR-InnerR)/(fN1-1);
   SetStepLength(stepLength,2*3.14159265/fN2);

   for(int i=n;i-->0;) fC1[i]=fC1[i]+InnerR;
   for (int i=0;i<n;i=i+fN1) {
      fV[i]=V1;
      fV[i+fN1-1]=V1;
      fIsFixed[i]=true;
      fIsFixed[i+fN1-1]=true;
   }
   SegmentId=SegmentId%Nseg;
   double segbegin=SegmentId*Z/Nseg;
   double segend=(SegmentId+1)*Z/Nseg;
   for (int i=fN1-1;i<n;i=i+fN1) {
     if(fC2[i]>=segbegin&&fC2[i]<=segend)
     {
        fV[i]=V0;
     }
     if(fC2[i]-segend<fdC2m[i]&&fC2[i]>segend)
     {
        fdC2m[i]=fC2[i]-segend;
     }
     if(-fC2[i]+segend<fdC2p[i]&&fC2[i]<segend)
     {
        fdC2p[i]=-fC2[i]+segend;
     }
   }
}

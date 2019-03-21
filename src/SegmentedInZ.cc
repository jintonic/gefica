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

   double stepLength=(OuterR-InnerR)/(N1-1);
   SetStepLength(stepLength,2*3.14159265/N2);

   for(int i=fN;i-->0;) C1[i]=C1[i]+InnerR;
   for (int i=0;i<fN;i=i+N1) {
      V[i]=Bias[1];
      V[i+N1-1]=Bias[1];
      fIsFixed[i]=true;
      fIsFixed[i+N1-1]=true;
   }
   SegmentId=SegmentId%Nseg;
   double segbegin=SegmentId*Z/Nseg;
   double segend=(SegmentId+1)*Z/Nseg;
   for (int i=N1-1;i<fN;i=i+N1) {
     if(C2[i]>=segbegin&&C2[i]<=segend)
     {
        V[i]=Bias[0];
     }
     if(C2[i]-segend<dC2m[i]&&C2[i]>segend)
     {
        dC2m[i]=C2[i]-segend;
     }
     if(-C2[i]+segend<dC2p[i]&&C2[i]<segend)
     {
        dC2p[i]=-C2[i]+segend;
     }
   }
}

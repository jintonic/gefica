#include "Units.h"
#include "SegmentedInPhi.h"
using namespace GeFiCa;

SegmentedInPhi::SegmentedInPhi(int nr, int np,
      const char *name, const char *title) : RhoPhi(nr, np, name, title),
   OuterR(2.5),InnerR(0.5),Nseg(6), SegmentId(1) {};
//______________________________________________________________________________
//
void SegmentedInPhi::InitializeGrid()
{
   if (InnerR>=OuterR) Fatal("InitializeGrid",
            "Inner R (%.1f) >= Outer R (%.1f)! Abort!", InnerR, OuterR);
   if (Nseg==0) Fatal("InitializeGrid",
            "Total number of segments cannot be zero! Abort!");
   
   if (Nseg<0)Nseg=-Nseg;
   if (SegmentId<0)SegmentId=-SegmentId;
   if(SegmentId>Nseg)SegmentId=SegmentId%Nseg;

   double steplength1=(OuterR-InnerR)/(N1-1);
   double steplength2=2*TMath::Pi()/(N2);
   double SegmentUpperBound=2*TMath::Pi()*SegmentId/Nseg;
   double SegmentLowerBound=2*TMath::Pi()*(SegmentId-1)/Nseg;
   SetStepLength(steplength1,steplength2);
   for(int i=0;i<fN;i++) {
      C1[i]+=InnerR;
      if(i%N1==0) {
         fIsFixed[i]=true;
         if(SegmentId==0)V[i]=Bias[1];
         else V[i]=Bias[0];
      } else if(i%N1==N1-1) { //need shift boundary as pointcontact
         fIsFixed[i]=true;
         if(SegmentId==0)V[i]=Bias[0];
         else if(C2[i]<=SegmentUpperBound&&C2[i]>=SegmentLowerBound) {
            V[i]=Bias[1];
         } else if(C2[i]>=SegmentUpperBound||C2[i]<=SegmentLowerBound) {
            V[i]=Bias[0];
         }
      } else {
         fIsFixed[i]=false;
         V[i]=0;
      }
   }
}

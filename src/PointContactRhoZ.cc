#include "Units.h"
#include "PointContactRhoZ.h"
using namespace GeFiCa;

void PointContactRhoZ::InitializeGrid()
{
   double stepLengthR=Radius/(N1-1);
   double stepLengthZ=Height/(N2-1);
   double RLowerBound=stepLengthR/2; // (2.0*Radius/(2*N1-1))/2;
   SetStepLength(stepLengthR,stepLengthZ);
   for(int i=fN;i-->0;) {
      C1[i]=C1[i]+RLowerBound;
      V[i]=(Bias[1]+Bias[0])/2;
      if(C1[i]>0&&C1[i]<PointContactR&&C2[i]<PointContactH) {
	    V[i]=Bias[0];
	    fIsFixed[i]=true;
      }
   }
   // set potential for electrodes
   for(int i=fN-1;i>=fN-N1;i--) {
      fIsFixed[i]=true;
      V[i]=Bias[1];
      if(C1[fN-1-i]>=0-0.001&&C1[fN-1-i]<=PointContactR+0.001) {
         V[fN-1-i]=Bias[0];
         fIsFixed[fN-1-i]=true;
      }
   }
   for(int i=0;i<fN-N1;i=i+N1) {
      fIsFixed[i+N1-1]=true;
      V[i+N1-1]=Bias[1];
   }

   for(int i=0;i<fN;i++) {
      //only change fdc1m when it is right close to bound
      if(C1[i]-PointContactR<dC1m[i]&&C1[i]>PointContactR
            &&C2[i]<PointContactH) dC1m[i]=C1[i]-PointContactR;
      //only change fdc2m when it is right close to bound
      if(C2[i]-PointContactH<dC2m[i]&&C2[i]>PointContactH
            &&C1[i]<PointContactR) dC2m[i]=C2[i]-PointContactH;
   }
}

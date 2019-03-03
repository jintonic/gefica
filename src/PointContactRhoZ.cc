#include "Units.h"
#include "PointContactRhoZ.h"
using namespace GeFiCa;

void PointContactRhoZ::Initialize()
{
   double stepLengthR=Radius/(n1-1);
   double stepLengthZ=Height/(n2-1);
   double RLowerBound=stepLengthR/2; // (2.0*Radius/(2*n1-1))/2;
   SetStepLength(stepLengthR,stepLengthZ);
   for(int i=n;i-->0;) {
      fC1[i]=fC1[i]+RLowerBound;
      fV[i]=(V1+V0)/2;
      if(fC1[i]>0&&fC1[i]<PointContactR&&fC2[i]<PointContactH) {
	    fV[i]=V0;
	    fIsFixed[i]=true;
      }
   }
   // set potential for electrodes
   for(int i=n-1;i>=n-n1;i--) {
      fIsFixed[i]=true;
      fV[i]=V1;
      if(fC1[n-1-i]>=0-0.001&&fC1[n-1-i]<=PointContactR+0.001) {
         fV[n-1-i]=V0;
         fIsFixed[n-1-i]=true;
      }
   }
   for(int i=0;i<n-n1;i=i+n1) {
      fIsFixed[i+n1-1]=true;
      fV[i+n1-1]=V1;
   }

   for(int i=0;i<n;i++) {
      //only change fdc1m when it is right close to bound
      if(fC1[i]-PointContactR<fdC1m[i]&&fC1[i]>PointContactR
            &&fC2[i]<PointContactH) fdC1m[i]=fC1[i]-PointContactR;
      //only change fdc2m when it is right close to bound
      if(fC2[i]-PointContactH<fdC2m[i]&&fC2[i]>PointContactH
            &&fC1[i]<PointContactR) fdC2m[i]=fC2[i]-PointContactH;
   }
}

#include "Units.h"
#include "PointContactRhoZ.h"
using namespace GeFiCa;

void PointContactRhoZ::InitializeGrid()
{
   double stepLengthR=Radius/(fN1-1);
   double stepLengthZ=Height/(fN2-1);
   double RLowerBound=stepLengthR/2; // (2.0*Radius/(2*fN1-1))/2;
   SetStepLength(stepLengthR,stepLengthZ);
   for(int i=fN;i-->0;) {
      fC1[i]=fC1[i]+RLowerBound;
      fV[i]=(V1+V0)/2;
      if(fC1[i]>0&&fC1[i]<PointContactR&&fC2[i]<PointContactH) {
	    fV[i]=V0;
	    fIsFixed[i]=true;
      }
   }
   // set potential for electrodes
   for(int i=fN-1;i>=fN-fN1;i--) {
      fIsFixed[i]=true;
      fV[i]=V1;
      if(fC1[fN-1-i]>=0-0.001&&fC1[fN-1-i]<=PointContactR+0.001) {
         fV[fN-1-i]=V0;
         fIsFixed[fN-1-i]=true;
      }
   }
   for(int i=0;i<fN-fN1;i=i+fN1) {
      fIsFixed[i+fN1-1]=true;
      fV[i+fN1-1]=V1;
   }

   for(int i=0;i<fN;i++) {
      //only change fdc1m when it is right close to bound
      if(fC1[i]-PointContactR<fdC1m[i]&&fC1[i]>PointContactR
            &&fC2[i]<PointContactH) fdC1m[i]=fC1[i]-PointContactR;
      //only change fdc2m when it is right close to bound
      if(fC2[i]-PointContactH<fdC2m[i]&&fC2[i]>PointContactH
            &&fC1[i]<PointContactR) fdC2m[i]=fC2[i]-PointContactH;
   }
}

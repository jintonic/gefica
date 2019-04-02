#include "Units.h"
#include "Detector.h"
using namespace GeFiCa;

Crystal::Crystal() : Height(1*cm),
   TopImpurity(1e10/cm3), BottomImpurity(1e10/cm3) {};
//______________________________________________________________________________
//
Detector::Detector(const char *name, const char *title)
   : Crystal(), TNamed(name,title) { Bias.push_back(0*volt); }
//______________________________________________________________________________
//
void Detector::CheckConfigurations()
{
   if (Height<=0) {
      Error("CheckConfigurations", "Height(%.1fcm)<=0!", Height/cm);
      abort();
   }
   if (Bias.size()<2) {
      Error("CheckConfigurations",
            "Number of electrodes == %zu!", Bias.size());
      abort();
   }
}

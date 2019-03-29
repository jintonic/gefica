#include "Units.h"
#include "TrueCoaxial.h"
using namespace GeFiCa;

TrueCoaxial::TrueCoaxial(const char *name, const char *title)
   : Detector(name, title), Radius(3*cm), BoreR(0.5*cm)
{ Bias.push_back(1*kV); }
//_____________________________________________________________________________
//
void TrueCoaxial::CheckConfigurations()
{
   Detector::CheckConfigurations();
   if (Radius<=0) {
      Error("CheckConfigurations", "Radius==%.1f!", Radius);
      abort();
   }
   if (BoreR<=0) {
      Error("CheckConfigurations", "BoreR==%.1f!", BoreR);
      abort();
   }
   if (BoreR>=Radius) {
      Error("CheckConfigurations", "BoreR(%.1f)>=Radius(%.1f)!",
            BoreR, Radius);
      abort();
   }
}

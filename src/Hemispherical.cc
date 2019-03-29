#include "Units.h"
#include "Hemispherical.h"
using namespace GeFiCa;

Hemispherical::Hemispherical(const char *name, const char *title)
   : Detector(name, title), PointContactR(2*mm), PointContactH(300*nm)
{ Bias.push_back(1*kV); }
//_____________________________________________________________________________
//
void Hemispherical::CheckConfigurations()
{
   Detector::CheckConfigurations();
   if (PointContactR<=0) {
      Error("CheckConfigurations", "PointContactR==%.1f!", PointContactR);
      abort();
   }
   if (PointContactH<=0) {
      Error("CheckConfigurations", "PointContactH==%.1f!", PointContactH);
      abort();
   }
   if (PointContactR>=Height) {
      Error("CheckConfigurations", "PointContactR(%.1f)>=Height(%.1f)!",
            PointContactR, Height);
      abort();
   }
   if (PointContactH>=Height) {
      Error("CheckConfigurations", "PointContactH(%.1f)>=Height(%.1f)!",
            PointContactH, Height);
      abort();
   }
}

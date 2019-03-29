#include "Units.h"
#include "Planar.h"
using namespace GeFiCa;

Planar::Planar(const char *name, const char *title)
   : Detector(name,title), Width(1*cm), Depth(1*cm)
{ Bias.push_back(1*kV); }
//_____________________________________________________________________________
//
void Planar::CheckConfigurations()
{
   Detector::CheckConfigurations();
   if (Width<=0) {
      Error("CheckConfigurations", "Width(%.1fcm)<=0!", Width/cm);
      abort();
   }
   if (Depth<=0) {
      Error("CheckConfigurations", "Depth(%.1fcm)<=0!", Depth/cm);
      abort();
   }
}

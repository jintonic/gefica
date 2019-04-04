#include "Units.h"
#include "Segmented.h"
using namespace GeFiCa;

Segmented::Segmented(const char *name, const char *title) :
   Detector(name, title), Radius(3.5*cm), BoreR(0.5*cm),
   Nphi(6), Nz(3), SegmentId(1)
{ Height=5*cm; Bias.push_back(2*kV); }
//______________________________________________________________________________
//
void Segmented::CheckConfigurations()
{
   if (Radius<=0) {
      Error("CheckConfigurations", "Radius==%.1f!", Radius);
      abort();
   }
   if (BoreR<=0) {
      Error("CheckConfigurations", "BoreR==%.1f!", BoreR);
      abort();
   }
   if (BoreR>=Radius) {
      Error("CheckConfigurations",
            "BoreR (%.1f) >= Radius (%.1f)!", BoreR, Radius);
      abort();
   }
   if (Nphi==0) {
      Error("CheckConfigurations",
         "Total number of segments in phi cannot be zero!");
      abort();
   }
   if (Nz==0) {
      Error("CheckConfigurations",
         "Total number of segments in z cannot be zero!");
      abort();
   }
   if (SegmentId>Nphi*Nz) {
      Error("CheckConfigurations",
            "SegmentId(%zu)>Nphi(%zu)*Nz(%zu)!",SegmentId, Nphi, Nz);
      abort();
   }
}

#include "Units.h"
#include "PointContact.h"
using namespace GeFiCa;

PointContact::PointContact(const char *name, const char *title)
   : Detector(name, title), Radius(3*cm),
   PointContactH(300*nm), PointContactR(1*mm),
   BoreH(0), BoreR(0), BoreTaperW(0), BoreTaperH(0),
   TaperW(1*mm), TaperH(1*mm), CornerW(1*mm), CornerH(1*mm),
   WrapAroundR(-1), GrooveW(0), GrooveH(0)
{ Height=5*cm; Bias.push_back(1*kV); }
//_____________________________________________________________________________
//
void PointContact::CheckConfigurations()
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
   if (BoreH<=0) {
      Error("CheckConfigurations", "BoreH==%.1f!", BoreH);
      abort();
   }
   if (BoreR<=0) {
      Error("CheckConfigurations", "BoreR==%.1f!", BoreR);
      abort();
   }
   if (BoreTaperW<=0) {
      Error("CheckConfigurations", "BoreTaperW==%.1f!", BoreTaperW);
      abort();
   }
   if (BoreTaperH<=0) {
      Error("CheckConfigurations", "BoreTaperH==%.1f!", BoreTaperH);
      abort();
   }
   if (TaperW<=0) {
      Error("CheckConfigurations", "TaperW==%.1f!", TaperW);
      abort();
   }
   if (TaperH<=0) {
      Error("CheckConfigurations", "TaperH==%.1f!", TaperH);
      abort();
   }
   if (CornerW<=0) {
      Error("CheckConfigurations", "CornerW==%.1f!", CornerW);
      abort();
   }
   if (CornerH<=0) {
      Error("CheckConfigurations", "CornerH==%.1f!", CornerH);
      abort();
   }
   if (GrooveW<=0) {
      Error("CheckConfigurations", "GrooveW==%.1f!", GrooveW);
      abort();
   }
   if (GrooveH<=0) {
      Error("CheckConfigurations", "GrooveH==%.1f!", GrooveH);
      abort();
   }
   if (GrooveH>=Height) {
      Error("CheckConfigurations", "GrooveH(%.1f)>=Height(%.1f)!",
            GrooveH, Height);
      abort();
   }
   if (WrapAroundR<0) WrapAroundR=Radius-TaperW;
   if (PointContactR+GrooveW>WrapAroundR) {
      Error("CheckConfigurations",
            "PointContactR(%.1f)+GrooveW(%.1f)>WrapAroundR(%.1f)!",
            PointContactR, GrooveH, WrapAroundR);
      abort();
   }
   if (PointContactH+BoreH>=Height) {
      Error("CheckConfigurations",
            "PointContactH(%.1f)+BoreH(%.1f)>=Height(%.1f)!",
            PointContactH, BoreH, Height);
      abort();
   }
   if (BoreR+BoreTaperW+CornerW>=Radius) {
      Error("CheckConfigurations",
            "BoreR(%.1f)+BoreTaperW(%.1f)+CornerW(%.1f)>=Radius(%.1f)!",
            BoreR, BoreTaperW, CornerW, Radius);
      abort();
   }
   if (TaperH+CornerH>=Height) {
      Error("CheckConfigurations", "TaperH(%.1f)+CornerH(%.1f)>=Height(%.1f)!",
            TaperH, CornerH, Height);
      abort();
   }
}

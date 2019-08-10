#include "Units.h"
#include "PointContact.h"
using namespace GeFiCa;

PointContact::PointContact(const char *name, const char *title) :
   Detector(name, title), Radius(3*cm),
   PointContactR(1*mm), PointContactH(300*nm),
   BoreH(0), BoreR(0), BoreTaperW(0), BoreTaperH(0),
   TaperW(1*mm), TaperH(1*mm), CornerW(1*mm), CornerH(1*mm),
   WrapAroundR(-1), GrooveW(0), GrooveH(0)
{ Height=5*cm; Bias.push_back(1*kV); }
//______________________________________________________________________________
//
void PointContact::CheckConfigurations()
{
   Detector::CheckConfigurations();
   if (Radius<=0) {
      Error("CheckConfigurations", "Radius==%.1f!", Radius);
      abort();
   }
   if (PointContactR<0) {
      Error("CheckConfigurations", "PointContactR==%.1f!", PointContactR);
      abort();
   }
   if (PointContactH<0) {
      Error("CheckConfigurations", "PointContactH==%.1f!", PointContactH);
      abort();
   }
   if (BoreH<0) {
      Error("CheckConfigurations", "BoreH==%.1f!", BoreH);
      abort();
   }
   if (BoreR<0) {
      Error("CheckConfigurations", "BoreR==%.1f!", BoreR);
      abort();
   }
   if (BoreTaperW<0) {
      Error("CheckConfigurations", "BoreTaperW==%.1f!", BoreTaperW);
      abort();
   }
   if (BoreTaperH<0) {
      Error("CheckConfigurations", "BoreTaperH==%.1f!", BoreTaperH);
      abort();
   }
   if (TaperW<0) {
      Error("CheckConfigurations", "TaperW==%.1f!", TaperW);
      abort();
   }
   if (TaperH<0) {
      Error("CheckConfigurations", "TaperH==%.1f!", TaperH);
      abort();
   }
   if (CornerW<0) {
      Error("CheckConfigurations", "CornerW==%.1f!", CornerW);
      abort();
   }
   if (CornerH<0) {
      Error("CheckConfigurations", "CornerH==%.1f!", CornerH);
      abort();
   }
   if (GrooveW<0) {
      Error("CheckConfigurations", "GrooveW==%.1f!", GrooveW);
      abort();
   }
   if (GrooveH<0) {
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
//______________________________________________________________________________
//
#include <TBox.h>
#include <TLine.h>
void PointContact::Draw(Option_t* option)
{
   TString opt(option); opt.ToLower();
   if (opt.Contains("half")) return; // not implemented yet
   // crystal outline
   TBox *out = new TBox(-Radius,0,Radius,Height); out->SetFillStyle(0);
   out->SetLineColor(kBlack); out->SetLineStyle(kDashed); out->Draw();
   // bore hole
   TBox *bore = new TBox(-BoreR,Height-BoreH,BoreR,Height);
   bore->SetLineColor(kBlack); bore->SetLineStyle(kDashed);
   bore->SetFillStyle(0); bore->Draw();
   // point contact
   TBox *pc = new TBox(-PointContactR,0,PointContactR,PointContactH);
   pc->SetLineColor(kBlack); pc->SetLineStyle(kDashed);
   pc->SetFillStyle(0); pc->Draw();
   // bottom tapers
   TLine *blt = new TLine(-Radius,TaperH,-Radius+TaperW,0);
   blt->SetLineColor(kBlack); blt->SetLineStyle(kDashed); blt->Draw();
   TLine *brt = new TLine(Radius-TaperW,0,Radius,TaperH);
   brt->SetLineColor(kBlack); brt->SetLineStyle(kDashed); brt->Draw();
   // top tapers
   TLine *tlt = new TLine(-Radius,Height-CornerH,-Radius+CornerW,Height);
   tlt->SetLineColor(kBlack); tlt->SetLineStyle(kDashed); tlt->Draw();
   TLine *trt = new TLine(Radius-CornerW,Height,Radius,Height-CornerH);
   trt->SetLineColor(kBlack); trt->SetLineStyle(kDashed); trt->Draw();
   // bore tapers
   TLine *lb = new TLine(-BoreTaperW-BoreR,Height,-BoreR,Height-BoreTaperH);
   lb->SetLineColor(kBlack); lb->SetLineStyle(kDashed); lb->Draw();
   TLine *rb = new TLine(BoreR,Height-BoreTaperH,BoreR+BoreTaperW,Height);
   rb->SetLineColor(kBlack); rb->SetLineStyle(kDashed); rb->Draw();
   // grove
   TBox *lg = new TBox(-WrapAroundR,0,-WrapAroundR+GrooveW,GrooveH);
   lg->SetLineColor(kBlack); lg->SetLineStyle(kDashed);
   lg->SetFillStyle(0); lg->Draw();
   TBox *rg = new TBox(WrapAroundR-GrooveW,0,WrapAroundR,GrooveH);
   rg->SetLineColor(kBlack); rg->SetLineStyle(kDashed);
   rg->SetFillStyle(0); rg->Draw();
}

#include "Units.h"
#include "SquarePointContact.h"
using namespace GeFiCa;

SquarePointContact::SquarePointContact(const char *name, const char *title) :
   Detector(name, title), Width(3*cm),Length(3*cm),
   PointContactW(1*mm),PointContactL(1*mm), PointContactH(300*nm),
   WrapAroundW(0),TaperW(0)
{ Height=5*cm; Bias.push_back(1*kV); }
//______________________________________________________________________________
//
void SquarePointContact::CheckConfigurations()
{
   //TODO
   Detector::CheckConfigurations();
   if (PointContactH<0) {
      Error("CheckConfigurations", "PointContactH==%.1f!", PointContactH);
      abort();
   }
}
//______________________________________________________________________________
//
#include <TBox.h>
#include <TLine.h>
void SquarePointContact::Draw(Option_t* option)
{
   TString opt(option); opt.ToLower();
   //TODO
   /*
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
   */
}

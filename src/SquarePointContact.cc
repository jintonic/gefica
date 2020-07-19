#include "Units.h"
#include "SquarePointContact.h"
using namespace GeFiCa;

SquarePointContact::SquarePointContact(const char *name, const char *title) :
   Detector(name, title), Width(1*cm),Length(1*cm),
   PointContactW(1*mm),PointContactL(1*mm), PointContactH(300*nm),
   WrapAroundW(2*mm), WrapAroundL(2*mm)
{ Height=1*cm; Bias.push_back(1*kV); }
//______________________________________________________________________________
//
void SquarePointContact::CheckConfigurations()
{
   Detector::CheckConfigurations();
   if (Width<=0) {
      Error("CheckConfigurations", "Width==%.1f!", Width); abort();
   }
   if (Length<=0) {
      Error("CheckConfigurations", "Length==%.1f!", Length); abort();
   }
   if (PointContactW<=0) {
      Error("CheckConfigurations", "PointContactW==%.1f!", PointContactW);
      abort();
   }
   if (PointContactL<=0) {
      Error("CheckConfigurations", "PointContactL==%.1f!", PointContactL);
      abort();
   }
   if (PointContactH<0) {
      Error("CheckConfigurations", "PointContactH==%.1f!", PointContactH);
      abort();
   }
   if (WrapAroundW<PointContactW || WrapAroundW>Width) {
      Error("CheckConfigurations", "WrapAroundW==%.1f!", WrapAroundW);
      Error("CheckConfigurations", "PointContactW==%.1f!", PointContactW);
      Error("CheckConfigurations", "Width==%.1f!", Width);
      abort();
   }
   if (WrapAroundL<PointContactL || WrapAroundL>Length) {
      Error("CheckConfigurations", "WrapAroundL==%.1f!", WrapAroundL);
      Error("CheckConfigurations", "PointContactL==%.1f!", PointContactL);
      Error("CheckConfigurations", "Length==%.1f!", Length);
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
   */
}

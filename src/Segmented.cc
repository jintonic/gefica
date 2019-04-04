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
   if (SegmentId==0) {
      Info("CheckConfigurations", "SegmentId==0, "
            "please use TrueCoaxial for the core electrode.");
      abort();
   }
}
//______________________________________________________________________________
//
#include <TLine.h>
#include <TEllipse.h>
void Segmented::Draw(Option_t* option)
{
   TString pointOfView(option); pointOfView.ToLower();
   if (pointOfView.Contains("top")) {
      double x=Radius*cos(2*Pi/Nphi), y=Radius*sin(2*Pi/Nphi);

      TLine *l1 = new TLine(-Radius,0,Radius,0);
      l1->SetLineColor(kBlack); l1->SetLineStyle(kDashed); l1->Draw();
      TLine *l2 = new TLine(-x,-y,x,y);
      l2->SetLineColor(kBlack); l2->SetLineStyle(kDashed); l2->Draw();
      TLine *l3 = new TLine(-x,y,x,-y);
      l3->SetLineColor(kBlack); l3->SetLineStyle(kDashed); l3->Draw();

      TEllipse *e1 = new TEllipse(0,0,BoreR,BoreR); e1->Draw();
      TEllipse *e2 = new TEllipse(0,0,Radius,Radius); e2->SetFillStyle(0);
      e2->Draw();
   }
}

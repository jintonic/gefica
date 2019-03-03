// compare numerical result to analytic calculation for a 1D true coaxial detector
using namespace GeFiCa;
void compare2analytic()
{
   // configure detector
   TrueCoaxial1D *num = new TrueCoaxial1D;
   num->V0=3000*volt; // core
   num->V1=0*volt; // outer surface electrode
   num->InnerR=0.5*cm;
   num->OuterR=3.0*cm;
   num->SetAverageImpurity(-3e9/cm3); // n-type

   // make a copy of the detector configuration
   TrueCoaxial1D *ana = (TrueCoaxial1D*) num->Clone("ana");

   // calculate potential using SOR method
   num->CalculatePotential(kSOR2);

   // fill grid with analytic result
   ana->CalculatePotential(kAnalytic);

   // prepare drawing style
   gROOT->SetStyle("Plain"); // pick up a good drawing style to modify
   gStyle->SetLegendBorderSize(0);
   gStyle->SetLegendFont(132);
   gStyle->SetLabelFont(132,"XY");
   gStyle->SetTitleFont(132,"XY");
   gStyle->SetLabelSize(0.05,"XY");
   gStyle->SetTitleSize(0.05,"XY");
   gStyle->SetPadRightMargin(0.01);
   gStyle->SetPadLeftMargin(0.12);
   gStyle->SetPadTopMargin(0.02);

   // generate graphics
   TTree *tn = num->GetTree();
   tn->Draw("v:c1");
   TGraph *gn = new TGraph(tn->GetSelectedRows(), tn->GetV2(), tn->GetV1());
   TTree *ta = ana->GetTree();
   ta->Draw("v:c1");
   TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());

   // compare numerical result to analytic calculation
   gn->SetMarkerColor(kBlue);
   gn->SetMarkerStyle(kCircle);
   gn->SetMarkerSize(0.8);
   ga->SetLineColor(kRed);
   gn->SetTitle(";Radius [cm];Potential [V]");
   gn->Draw("ap");
   ga->Draw("l");

   TLegend *l = new TLegend(0.5,0.65,0.8,0.8);
   l->AddEntry(ga,"Analytic","l");
   l->AddEntry(gn,"SOR2","p");
   l->Draw();

   gPad->Print("tc1.png");
}

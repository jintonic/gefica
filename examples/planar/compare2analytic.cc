// compare numerical result to analytic calculation for a 1D planar detector
using namespace GeFiCa;
void compare2analytic()
{
   // configure detector
   Planar1D *num = new Planar1D;
   num->V0=0*volt;
   num->V1=800*volt;
   num->Thickness=1*cm;
   num->SetAverageImpurity(1e10/cm3);

   // make a copy of the detector configuration
   Planar1D *ana = (Planar1D*) num->Clone("ana");

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
   gn->SetTitle(";Thickness [cm];Potential [V]");
   gn->GetXaxis()->SetRangeUser(0,1);
   gn->GetYaxis()->SetRangeUser(0,800);
   //gn->GetYaxis()->SetTitleOffset(1.2);
   gn->Draw("ap");

   ga->SetLineColor(kRed);
   ga->Draw("l");

   TLegend *l = new TLegend(0.2,0.6,0.4,0.8);
   l->AddEntry(ga,"Analytic","l");
   l->AddEntry(gn,"SOR2","p");
   l->Draw();

   gPad->Print("planar1d.png");
}

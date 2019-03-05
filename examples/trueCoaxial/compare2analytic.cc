using namespace GeFiCa;
// compare numerical result to analytic calculation for a 1D true coaxial detector
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

   // generate graphics
   gStyle->SetPadRightMargin(0.02);

   TTree *tn = num->GetTree();
   tn->Draw("v:c1","","goff");
   TGraph *gn = new TGraph(tn->GetSelectedRows(), tn->GetV2(), tn->GetV1());
   TTree *ta = ana->GetTree();
   ta->Draw("v:c1","","goff");
   TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());

   // compare numerical result to analytic calculation
   gn->SetMarkerColor(kBlue);
   gn->SetMarkerStyle(kCircle);
   gn->SetMarkerSize(0.8);
   ga->SetLineColor(kRed);
   gn->SetTitle(";Radial position [cm];Potential [V]");
   gn->Draw("ap");
   ga->Draw("l");

   TLegend *l = new TLegend(0.5,0.65,0.8,0.8);
   l->AddEntry(ga,"Analytic","l");
   l->AddEntry(gn,"SOR","p");
   l->Draw();

   gPad->Print("tc1.png");
}

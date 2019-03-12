using namespace GeFiCa;
//Compare numerical result to analytic calculation for 1D spherical detector
void compare2analytic()
{
   // configure detector
   Sphere1D *num=new Sphere1D;
   num->InnerR=0.3*cm;
   num->OuterR=1*cm;
   num->SetAverageImpurity(3e9/cm3);
   num->V0=900*volt;
   num->V1=0*volt;
   num->Dump();
   cout<<"press any key to continue"<<endl; cin.get();

   // make a copy of the detector configuration
   Sphere1D *ana = (Sphere1D*) num->Clone("ana");

   // calculate potential using SOR method
   num->CalculatePotential();

   // fill grid with analytic result
   ana->FillGridWithAnalyticResult();

   // generate graphics
   TTree *tn = num->GetTree();
   tn->Draw("v:c1","","goff");
   TGraph *gn = new TGraph(tn->GetSelectedRows(), tn->GetV2(), tn->GetV1());
   TTree *ta = ana->GetTree();
   ta->Draw("v:c1","","goff");
   TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());

   // compare numerical result to analytic calculation
   gStyle->SetPadLeftMargin(0.11); gStyle->SetPadRightMargin(0.01);
   gn->SetMarkerColor(kBlue);
   gn->SetMarkerStyle(kCircle);
   gn->SetMarkerSize(0.8);
   gn->SetTitle(";Radial position [cm];Potential [V]");
   gn->GetXaxis()->SetRangeUser(0,1);
   gn->GetYaxis()->SetRangeUser(0,900);
   gn->Draw("ap");

   ga->SetLineColor(kRed);
   ga->Draw("l");

   TLegend *l = new TLegend(0.6,0.6,0.8,0.8);
   l->AddEntry(ga,"Analytic","l");
   l->AddEntry(gn,"SOR","p");
   l->Draw();

   gPad->Print("s1d.png");
}

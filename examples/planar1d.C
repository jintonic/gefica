// compare numerical result to analytic calculation for 1D planar detector
using namespace GeFiCa;
void planar1d()
{
   // define detector
   Planar1D *detector = new Planar1D(101);
   detector->MaxIterations=1e5;
   detector->Csor=1.95;
   detector->UpperBound=1*cm;
   detector->V0=0*volt;
   detector->V1=800*volt;
   detector->SetAverageImpurity(1e10/cm3);
   detector->Dump();
   cout<<"press any key to continue"<<endl;
   cin.get();

   // calculate fields with two different methods
   detector->CalculatePotential(kSOR2);
   detector->SaveField("planar1dSOR2.root");
   detector->CalculatePotential(kAnalytic);
   detector->SaveField("planar1dTRUE.root");

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

   // compare numerical result to analytic calculation
   TChain *tn = new TChain("t");
   tn->Add("planar1dSOR2.root");
   tn->Draw("v:c1");
   TGraph *gn = new TGraph(tn->GetSelectedRows(), tn->GetV2(), tn->GetV1());

   TChain *ta = new TChain("t");
   ta->Add("planar1dTRUE.root");
   ta->Draw("v:c1");
   TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());

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

   // calculate capacitance
   double c = detector->GetCapacitance()/pF;
   cout<<"capacitance is "<<c<<" pF per cm2"<<endl;
}

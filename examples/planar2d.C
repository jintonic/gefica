// compare numerical result to analytic calculation for 2D planar detector
{
   GeFiCa::Planar2D *detector2 = new GeFiCa::Planar2D(101,101);
   detector2->MaxIterations=1e5;
   detector2->Csor=1.995;
   detector2->V0=800*GeFiCa::volt;
   detector2->V1=0*GeFiCa::volt;
   detector2->SetAverageImpurity(1e10/GeFiCa::cm3);
   detector2->Dump();
   cout<<"press any key to continue"<<endl;
   cin.get();

   //use analytic method from 1D planar
   GeFiCa::Planar1D *detector = new GeFiCa::Planar1D(101);
   detector->V0=800*GeFiCa::volt;
   detector->V1=0*GeFiCa::volt;
   detector->SetAverageImpurity(1e10/GeFiCa::cm3);

   // calculate fields with two different methods
   detector2->CalculatePotential(GeFiCa::kSOR2);
   detector2->SaveField("planar2dSOR2.root");
   detector->CalculatePotential(GeFiCa::kAnalytic);
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

   // generate graphics
   TChain *tn = new TChain("t");
   tn->Add("planar2dSOR2.root");
   tn->Draw("v:c1");
   TGraph *gn = new TGraph(tn->GetSelectedRows(), tn->GetV2(), tn->GetV1());

   TChain *ta = new TChain("t");
   ta->Add("planar1dTRUE.root");
   ta->Draw("v:c1");
   TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());

   // make final plot
   gn->SetMarkerColor(kBlue);
   gn->SetMarkerStyle(kCircle);
   gn->SetMarkerSize(0.8);
   ga->SetLineColor(kRed);
   gn->SetTitle(";Thickness [cm];Potential [V]");
   gn->Draw("ap");
   ga->Draw("l");

   TLegend *l = new TLegend(0.7,0.6,0.9,0.8);
   l->AddEntry(ga,"Analytic","l");
   l->AddEntry(gn,"SOR2","p");
   l->Draw();

   gPad->Print("planar2d.png");

   // calculate capacitance
   double c = detector->GetCapacitance()/GeFiCa::pF;
   cout<<"capacitance is "<<c<<" pF per cm2"<<endl;
}

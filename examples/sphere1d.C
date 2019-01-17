// compare numerical result to analytic calculation for 1D Sphere detector
{
   GeFiCa::Sphere1D *detector=new GeFiCa::Sphere1D(101);
   detector->MaxIterations=1e5;
   detector->Csor=1.9;
   detector->V1=0*GeFiCa::volt;
   detector->V0=900*GeFiCa::volt;
   detector->InnerRadius=0.3*GeFiCa::cm;
   detector->OuterRadius=1*GeFiCa::cm;
   detector->SetAverageImpurity(1e10/GeFiCa::cm3);
   detector->Dump();
   cout<<"press any key to continue"<<endl;
   cin.get();

   // calculate fields with two different methods
   detector->CalculatePotential(GeFiCa::kSOR2);
   detector->SaveField("sphere1dSOR2.root");
   detector->CalculatePotential(GeFiCa::kAnalytic);
   detector->SaveField("sphere1dTRUE.root");

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
   tn->Add("sphere1dSOR2.root");
   tn->Draw("p:c1");
   TGraph *gn = new TGraph(tn->GetSelectedRows(), tn->GetV2(), tn->GetV1());

   TChain *ta = new TChain("t");
   ta->Add("sphere1dTRUE.root");
   ta->Draw("p:c1");
   TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());

   gn->SetMarkerColor(kBlue);
   gn->SetMarkerStyle(kCircle);
   gn->SetMarkerSize(0.8);
   gn->SetTitle(";Radius [cm];Potential [V]");
   gn->GetXaxis()->SetRangeUser(0,3);
   gn->GetYaxis()->SetRangeUser(0,900);
   gn->Draw("ap");

   ga->SetLineColor(kRed);
   ga->Draw("l");

   TLegend *l = new TLegend(0.6,0.6,0.8,0.8);
   l->AddEntry(ga,"Analytic","l");
   l->AddEntry(gn,"SOR2","p");
   l->Draw();

   gPad->Print("sphere1d.png");
}

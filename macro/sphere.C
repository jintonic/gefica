{
  GeFiCa::Sphere1D *detector=new GeFiCa::Sphere1D(101);
  detector->MaxIterations=1e5;
  detector->Csor=1.95;
  detector->Vpos=2000*GeFiCa::volt;
  detector->Vneg=0*GeFiCa::volt;
  detector->SetImpurity(1e10/GeFiCa::cm3);
  detector->CalculateField(GeFiCa::kSOR2);
  detector->SaveField("sphere1dSOR2.root");
  detector->CalculateField(GeFiCa::kAnalytic);
  detector->SaveField("sphere1dTrue.root");

   // generate graphics
   TChain *tn = new TChain("t");
   tn->Add("sphere1dSOR2.root");
   tn->Draw("p:c1");
   TGraph *gn = new TGraph(tn->GetSelectedRows(), tn->GetV2(), tn->GetV1());

   TChain *ta = new TChain("t");
   ta->Add("sphere1dTrue.root");
   ta->Draw("p:c1");
   TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());

   // make final plot
   gn->SetMarkerColor(kBlue);
   gn->SetMarkerStyle(6);
   gn->SetMarkerSize(0.8);
   ga->SetLineColor(kRed);
   gn->SetTitle(";Thickness [cm];Potential [V]");
   gn->Draw("ap");
   ga->Draw("l");

}

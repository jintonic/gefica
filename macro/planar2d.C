{
   GeFiCa::Planar2D *detector2 = new GeFiCa::Planar2D(101,101);
   detector2->MaxIterations=1e5;
   detector2->Csor=1.995;
   detector2->V0=2000*GeFiCa::volt;
   detector2->V1=0*GeFiCa::volt;
   detector2->Impurity="1e10";
   detector2->CalculateField(GeFiCa::kSOR2);
   detector2->SaveField("planar2dSOR2.root");

   // calculate fields
   GeFiCa::Planar1D *detector = new GeFiCa::Planar1D(101);
   detector->V0=2000*GeFiCa::volt;
   detector->V1=0*GeFiCa::volt;
   detector->SetImpurity(1e10/GeFiCa::cm3);
   detector->CalculateField(GeFiCa::kAnalytic);
   detector->SaveField("planar1dTrue.root");

   // generate graphics
   TChain *tn = new TChain("t");
   tn->Add("planar2dSOR2.root");
   tn->Draw("p:c1");
   TGraph *gn = new TGraph(tn->GetSelectedRows(), tn->GetV2(), tn->GetV1());

   TChain *ta = new TChain("t");
   ta->Add("planar1dTrue.root");
   ta->Draw("p:c1");
   TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());

   // make final plot
   gn->SetMarkerColor(kBlue);
   gn->SetMarkerStyle(8);
   gn->SetMarkerSize(0.8);
   ga->SetLineColor(kRed);
   gn->SetTitle(";Thickness [cm];Potential [V]");
   gn->Draw("ap");
   ga->Draw("l");
}

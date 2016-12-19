{
   GeFiCa::TrueCoaxial1D *detector = new GeFiCa::TrueCoaxial1D(101);
   detector->Vpos=0*GeFiCa::volt;
   detector->Vneg=3000*GeFiCa::volt;
   detector->MaxIterations=1e5;
   detector->SetImpurity(1e10/GeFiCa::cm3);
   detector->Csor=1.95;
   d=detector->CalculateField(GeFiCa::kSOR2);
   detector->SaveField("trueCoaxial1d.root");
   detector->CalculateField(GeFiCa::kAnalytic);
   detector->SaveField("trueCoaxial1dTrue.root");

   // generate graphics
   TChain *tn = new TChain("t");
   tn->Add("trueCoaxial1dTrue.root");
   tn->Draw("p:c1");
   TGraph *gn = new TGraph(tn->GetSelectedRows(), tn->GetV2(), tn->GetV1());

   TChain *ta = new TChain("t");
   ta->Add("trueCoaxial1d.root");
   ta->Draw("p:c1");
   TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());

   // make final plot
   gn->SetMarkerColor(kBlue);
   gn->SetMarkerStyle(6);
   gn->SetMarkerSize(0.8);
   ga->SetLineColor(kRed);
   gn->SetTitle(";Thickness [cm];Potential [V]");
   ga->Draw("al");
   gn->Draw("p");

}

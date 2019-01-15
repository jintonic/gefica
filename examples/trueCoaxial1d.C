{
   GeFiCa::TrueCoaxial1D *detector = new GeFiCa::TrueCoaxial1D(332);
   detector->V0=0*GeFiCa::volt;
   detector->V1=5500*GeFiCa::volt;
   detector->InnerRadius=0.15;
   detector->OuterRadius=3.05;

   detector->MaxIterations=1e5;
   detector->SetAverageImpurity(-0.318e10/GeFiCa::cm3);
   detector->Csor=1.95;
   detector->CalculatePotential(GeFiCa::kSOR2);
   detector->SaveField("trueCoaxial1d.root");
   detector->CalculatePotential(GeFiCa::kAnalytic);
   detector->SaveField("trueCoaxial1dTrue.root");

   // generate graphics
   TChain *tn = new TChain("t");
   tn->Add("trueCoaxial1d.root");
   tn->Draw("p:c1*10");
   TGraph *gn = new TGraph(tn->GetSelectedRows(), tn->GetV2(), tn->GetV1());


   TChain *ta = new TChain("t");
   ta->Add("trueCoaxial1dTrue.root");
   ta->Draw("p:c1*10");
   TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());

   // make final plot
   gn->SetMarkerColor(kBlue);
   gn->SetMarkerStyle(kCircle);
   gn->SetMarkerSize(0.8);
   ga->SetLineColor(kRed);
   gn->SetTitle(";Thickness [cm];Potential [V]");
   gn->Draw("ap");
   ga->Draw("l");

}

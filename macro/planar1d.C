{
   GeFiCa::Planar1D *detector = new GeFiCa::Planar1D(101);
   detector->MaxIterations=1e5;
   detector->Csor=1.9;
   detector->SetImpurity(1e10/GeFiCa::cm3);
   detector->CalculateField(GeFiCa::kSOR2);
   detector->SaveField("planar1dSOR2.root");
   detector->CalculateField(GeFiCa::kAnalytic);
   detector->SaveField("planar1dTrue.root");

   TChain *t1 = new TChain("t");
   t1->Add("planar1dSOR2.root");
   t1->Draw("p:c1");

   TChain *t2 = new TChain("t");
   t2->Add("planar1dTrue.root");
   t2->Draw("p:c1","","same");
}

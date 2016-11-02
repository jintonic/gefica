{
   GeFiCa::Planar2D *detector = new GeFiCa::Planar2D(101,101);
   detector->MaxIterations=1e5;
   detector->Csor=1.9;
   detector->SetImpurity(1e10/GeFiCa::cm3);
   detector->CalculateField(GeFiCa::kSOR2);
   detector->SaveField("planar2d.root");

   TChain *t = new TChain("t");
   t->Add("planar2d.root");
   t->Draw("p:c1:c2");
}

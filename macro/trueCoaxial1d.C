{
   GeFiCa::TrueCoaxial1D *detector = new GeFiCa::TrueCoaxial1D(1001);
   detector->MaxIterations=1e5;
   detector->Csor=1.9;
   detector->CalculateField(GeFiCa::kSOR2);
   detector->SaveField("trueCoaxial1d.root");
   detector->CalculateField(GeFiCa::kAnalytic);
   detector->SaveField("trueCoaxial1dTrue.root");

   TChain *t1=new TChain("t");
   t1->Add ("trueCoaxial1d.root");
   TChain *t2=new TChain("t");
   t2->Add("trueCoaxial1dTrue.root");
   t1->Draw("p:c1");
   t2->Draw("p:c1","","same");
}

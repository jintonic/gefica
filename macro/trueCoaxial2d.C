{
   GeFiCa::TrueCoaxial2D *detector = new GeFiCa::TrueCoaxial2D(101,101);
   detector->MaxIterations=1e5;
   detector->Csor=1.95;
   detector->CalculateField(GeFiCa::kSOR2);
   detector->SaveField("trueCoaxial2d.root");

   TChain *t=new TChain("t");
   t->Add("trueCoaxial2d.root");
   t->Draw("p:c1:c2");
}

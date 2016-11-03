{
   GeFiCa::PointContactXY *detector = new GeFiCa::PointContactXY(1001,1001);
   detector->MaxIterations=1e5;
   detector->Csor=1.9;
   detector->SetImpurity(1e10/GeFiCa::cm3);
   detector->YUpperBound=5;
   detector->XUpperBound=13;
   detector->PointBegin=5.5;
   detector->PointEnd=6.5;
   detector->cathode_voltage=1000;
   detector->Initialize();

   detector->CalculateField(GeFiCa::kSOR2);
   detector->SaveField("pointContact2d.root");

   TChain *t = new TChain("t");
   t->Add("pointContact2d.root");
   t->Draw("c2:c1:p", "", "colz");
}

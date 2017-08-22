{
   // calculate fields
   GeFiCa::Planar1D *detector = new GeFiCa::Planar1D(101);
   detector->MaxIterations=1e5;
   detector->Csor=1.95;
   detector->V1=2000*GeFiCa::volt;
   detector->V0=0*GeFiCa::volt;
   detector->SetImpurity(1e10/GeFiCa::cm3);
   detector->CalculateField(GeFiCa::kSOR2);
   detector->SaveField("planar1dSOR2.root");
   detector->CalculateField(GeFiCa::kAnalytic);
   //save to local .root file
   detector->SaveField("planar1dTrue.root");
}

using namespace GEFICA;

void polar1d()
{
   Polar1d *detector = new Polar1d(11);
   detector->Thickness=10*cm;
   detector->MaxIterations=1e5;
   detector->Csor=1;
   detector->Create(0.3,3);
   detector->SetVoltage(0*volt, 2000*volt);
   detector->SetImpurity(1e10/cm3);
   detector->CalculateField(EMethod::kRK2);
   detector->SaveField("planar1dRK2.root");
   detector->CalculateField(EMethod::kAnalytic);
   detector->SaveField("planar1dTrue.root");
}

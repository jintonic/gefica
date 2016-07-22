using namespace GEFICA;

void planar1d()
{
   PlanarX *detector = new PlanarX(5001);
   detector->Thickness=10*cm;
   detector->MaxIterations=1e6;
   detector->Csor=1.9;
   detector->SetVoltage(0*volt, 2000*volt);
   detector->SetImpurity(1e10/cm3);
   detector->CalculateField(EMethod::kRK2);
   detector->SaveField("planar1dRK2.root");
   detector->CalculateField(EMethod::kAnalytic);
   detector->SaveField("planar1dTrue.root");
}

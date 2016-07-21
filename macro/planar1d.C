using namespace GEFICA;

void planar1d()
{
   PlanarX *detector = new PlanarX(500);
   detector->Thickness=10*cm;
   detector->MaxIterations=1e7;
   detector->SetVoltage(0*volt, 2000*volt);
   detector->SetImpurity(1e10/cm3);
   detector->CalculateField(EMethod::kRK4);
   detector->SaveField("planar1dRK2.root");
   detector->CalculateField(EMethod::kAnalytic);
   detector->SaveField("planar1dTrue.root");
}

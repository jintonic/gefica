using namespace GEFICA;

void planar1d()
{
   PlanarX *detector = new PlanarX;
   detector->Thickness=10*cm;
   detector->SetVoltage(0*volt, 2000*volt);
   detector->SetImpurity(1e10/cm3);
   detector->CalculateField(EMethod::kRK2);
   detector->SaveField("planar1dRK2.root");
   detector->CalculateField(EMethod::kAnalytic);
   detector->SaveField("planar1dTrue.root");
}

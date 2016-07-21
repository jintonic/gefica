using namespace GEFICA;

void planar1d(const char* file="planar1d.root")
{
   PlanarX *detector = new PlanarX;
   detector->Thickness=10*cm;
   detector->SetVoltage(0, 2000*volt);
   detector->SetImpurity(1e10/cm3);
   detector->CalculateField(EMethod::kRK2);
   detector->SaveField(file);
}

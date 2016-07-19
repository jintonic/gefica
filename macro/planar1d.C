#include <X.h>
using namespace GEFICA;

void planar1d(const char* fin="input.root")
{
   PlanarX *field = new PlanarX(thickness=10/*cm*/);
   field->n=101;
   field->SetVoltage(2000/*volt*/);
   field->Calculate(MaxIterations=10000);
   field->Save("output.root");
}

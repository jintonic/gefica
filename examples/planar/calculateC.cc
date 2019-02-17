// calculate capacitance in unit area of an ideal 1D planar detector
using namespace GeFiCa;
void calculateC()
{
   // configure detector
   Planar1D *detector = new Planar1D;
   detector->Thickness=1*cm;
   detector->V0=0*volt;
   detector->V1=800*volt;

   // calculate capacitance
   double c = detector->GetCapacitance()/pF;
   cout<<"capacitance is "<<c<<" pF/cm2"<<endl;
}

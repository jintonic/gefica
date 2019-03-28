// calculate capacitance in unit area of an ideal planar detector
using namespace GeFiCa;
void calculateC()
{
   Planar detector;
   detector.Height=1*cm;
   detector.Bias[0]=0*volt;
   detector.Bias[1]=800*volt;
   detector.SetAverageImpurity(1e10/cm3);

   X grid;
   grid.GetBoundaryConditionFrom(detector);
   double c = grid.GetC()/pF;
   cout<<"capacitance is "<<c<<" pF/cm2"<<endl;
}

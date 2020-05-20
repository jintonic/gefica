using namespace GeFiCa;
// calculate and save fields of an Inverted Coaxial Point-Contact (ICPC) HPGe
void calculateFields(const char *output="ICPC.root")
{
   SquarePointContact detector;
   detector.Bias[0]=100*kV; // bias on point contact
   detector.Bias[1]=0*kV; // ground outer contact

   detector.Width=1.45*cm; detector.Height=1.05*cm; detector.Length=1*cm;

   detector.PointContactW=1*mm; detector.PointContactH=1*mm;
   detector.PointContactL=5*mm;

   detector.TaperW=0.1*cm;
   detector.WrapAroundW=0.5*cm;

   detector.BottomImpurity=3e9/cm3; detector.TopImpurity=7e9/cm3;

   int nx=50; 
   int ny=50; // precision: 0.1 mm
   int nz=50; // precision: 0.1 mm
   XYZ grid(nx,ny,nz);
   grid.SetupWith(detector);
   grid.RelaxationFactor=1.94;
   grid.SuccessiveOverRelax();
   
   TFile file(output,"recreate");
   detector.Write();
   grid.Write();
   file.Close();
}

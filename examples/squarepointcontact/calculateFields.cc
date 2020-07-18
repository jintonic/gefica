using namespace GeFiCa;
// calculate and save fields of an Inverted Coaxial Point-Contact (ICPC) HPGe
void calculateFields(const char *output="SPC.root")
{
   SquarePointContact detector;
   detector.Bias[0]=70; // bias on point contact
   detector.Bias[1]=0; // ground outer contact

   detector.Width=1.8*cm; detector.Height=1.0*cm; detector.Length=1.8*cm;

   detector.PointContactW=0.6*mm; detector.PointContactH=0.1*mm;
   detector.PointContactL=0.6*mm;

   detector.TaperW=0.3*cm;
   //detector.WrapAroundW=0.3*cm;

   detector.BottomImpurity=4e9/cm3; detector.TopImpurity=4e9/cm3;

   int nx=50; 
   int ny=50; // precision: 0.1 mm
   int nz=50; // precision: 0.1 mm
   XYZ grid(nx,ny,nz);
   grid.SetupWith(detector);
   grid.RelaxationFactor=1.84;
   grid.SuccessiveOverRelax();
   
   TFile file(output,"recreate");
   detector.Write();
   grid.Write();
   file.Close();
}

using namespace GeFiCa;
// calculate and save fields of an Inverted Coaxial Point-Contact (ICPC) HPGe
void calculateFields(const char *output="ICPC.root")
{
   PointContact detector;
   detector.Bias[0]=-2.5*kV; // bias on point contact
   detector.Bias[1]=0; // ground outer contact

   detector.Radius=3.45*cm; detector.Height=5.05*cm;

   detector.PointContactR=1.4*mm; detector.PointContactH=0.1*mm;

   detector.TaperW=3*mm; detector.TaperH=3*mm;
   detector.CornerW=3*mm; detector.CornerH=3*mm;

   detector.BoreH=2*cm; detector.BoreR=8*mm;
   detector.BoreTaperH=3*mm; detector.BoreTaperW=3*mm;

   detector.WrapAroundR=2.5*cm;
   detector.GrooveH=3*mm; detector.GrooveW=3*mm;

   detector.BottomImpurity=3e9/cm3; detector.TopImpurity=7e9/cm3;

   int nr=345*2; // no point on r=0, please
   int nz=505+1; // precision: 0.1 mm
   RhoZ grid(nr,nz);
   grid.SetupWith(detector);
   grid.RelaxationFactor=1.994;
   grid.SuccessiveOverRelax();
   
   TFile *file = new TFile(output,"recreate");
   detector.Write();
   grid.Write();
   file->Close();
}

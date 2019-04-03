using namespace GeFiCa;
// calculate and save fields of a PPC
void calculateFields(const char *output="ppc.root")
{
   PointContact detector;
   detector.Bias[0]=-2.5*kV; // bias on point contact
   detector.Bias[1]=0; // ground outer contact

   detector.Radius=3.45*cm; detector.Height=5.05*cm;

   detector.PointContactR=1.40*mm; detector.PointContactH=2.10*mm;

   detector.WrapAroundR=1.20*cm;
   detector.TaperW=0.1*cm; detector.TaperH=0.1*cm;
   detector.CornerW=0.1*cm; detector.CornerH=0.1*cm;

   detector.BoreH=2*cm; detector.BoreR=2*cm;
   detector.BoreTaperH=2*mm; detector.BoreTaperW=2*mm;

   detector.GrooveH=1*mm; detector.GrooveW=1*mm;

   detector.BottomImpurity=-3.18e9/cm3; detector.TopImpurity=-1.9175e9/cm3;

   RhoZ grid;
   grid.GetBoundaryConditionFrom(detector);
   grid.RelaxationFactor=1.995;
   grid.SuccessiveOverRelax();
   
   TFile *file = new TFile(output,"recreate");
   detector.Write();
   grid.GetTree()->Write();
   file->Close();
}

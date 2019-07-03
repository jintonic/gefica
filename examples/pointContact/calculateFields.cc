using namespace GeFiCa;
// calculate and save fields of a PPC
void calculateFields(const char *output="ppc.root")
{
   PointContact detector;
   detector.Bias[0]=0*kV; // bias on point contact
   detector.Bias[1]=2500; // ground outer contact

   detector.Radius=3.45*cm; detector.Height=5.05*cm;

//   detector.PointContactR=1.40*mm; detector.PointContactH=2.10*mm;
//   detector.WrapAroundR=1.50*mm;
//   detector.TaperW=0.4*cm; detector.TaperH=0.4*cm;
//   detector.CornerW=0.4*cm; detector.CornerH=0.4*cm;

//   detector.BoreH=1*cm; detector.BoreR=1*cm;
//   detector.BoreTaperH=3*mm; detector.BoreTaperW=3*mm;
//   detector.GrooveH=0.001*mm; detector.GrooveW=.5*mm;

   detector.BottomImpurity=0/cm3; detector.TopImpurity=0/cm3;

   RhoZ grid(500,500);
   grid.GetBoundaryConditionFrom(detector);
//   grid.MaxIterations=2; ///< maximal iterations of SOR to be performed
   grid.RelaxationFactor=1.995;
   grid.SuccessiveOverRelax();
   
   TFile *file = new TFile(output,"recreate");
   detector.Write();
   grid.Write();
   file->Close();
}

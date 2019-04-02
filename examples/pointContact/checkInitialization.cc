// show default grid setup of 2D point contact detectors
{
   GeFiCa::PointContact detector;
   detector.Dump(); // call TObject::Dump() to print data members
   detector.BoreR=5*GeFiCa::mm; detector.BoreH=3*GeFiCa::cm; // add bore hole
   detector.BoreTaperW=2*GeFiCa::mm; detector.BoreTaperH=2*GeFiCa::mm; // add taper
   detector.WrapAroundR=2*GeFiCa::cm; // add wrap arround
   detector.GrooveW=2*GeFiCa::mm; detector.GrooveH=2*GeFiCa::mm; // add groove

   GeFiCa::RhoZ grid;
   grid.GetBoundaryConditionFrom(detector);
   TTree *t = grid.GetTree(); // create a ROOT tree for quick check
   gStyle->SetPadRightMargin(0.12); // give enough space for color palette
   t->Draw("c2:c1:v","","colz"); // visualize initial setup
}

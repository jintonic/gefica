// show default grid setup of 2D point contact detectors
{
   GeFiCa::PointContact detector;
   detector.BoreR=5*GeFiCa::mm; detector.BoreH=3*GeFiCa::cm; // add bore hole
   detector.BoreTaperW=2*GeFiCa::mm; detector.BoreTaperH=2*GeFiCa::mm;
   detector.TaperW=5*GeFiCa::mm; detector.TaperH=4*GeFiCa::mm;
   detector.CornerW=5*GeFiCa::mm; detector.CornerH=4*GeFiCa::mm;
   detector.WrapAroundR=2*GeFiCa::cm; // add wrap arround
   detector.GrooveW=2*GeFiCa::mm; detector.GrooveH=2*GeFiCa::mm;
   detector.Dump(); // call TObject::Dump() to print data members

   GeFiCa::RhoZ grid(40,21);
   grid.GetBoundaryConditionFrom(detector);
   TTree *t = grid.GetTree(); // create a ROOT tree for quick check
   gStyle->SetPadRightMargin(0.12); // give enough space for color palette
   t->Draw("c2:c1:p1","","colz"); // visualize initial setup
   detector.Draw();
   TCanvas c1;
   t->Draw("c2:c1:m1","","colz"); // visualize initial setup
   detector.Draw();
}

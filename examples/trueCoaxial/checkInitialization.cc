// show default configuration & grid setup of ideal true coaxial detectors
{
   GeFiCa::TrueCoaxial detector; // construct an ideal true coaxial detector
   detector.Dump(); // call TObject::Dump() to print data members
   GeFiCa::Rho grid;
   grid.GetBoundaryConditionFrom(detector);
   TTree *t = grid.GetTree(); // create a ROOT tree for quick check
   t->Draw("v:c1");
}

// show default configuration & grid setup of 1d planar detectors
{
   GeFiCa::Planar detector;
   detector.Dump(); // call TObject::Dump() to print default detector setup

   GeFiCa::X grid;
   detector.Configure(grid);
   
   TTree *t = grid.GetTree(); // create a ROOT tree for quick investigation
   t->Draw("v:c1");
}

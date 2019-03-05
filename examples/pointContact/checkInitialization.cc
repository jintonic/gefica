// show default grid setup of 2D point contact detectors
{
   GeFiCa::PointContactDZ pcdz;
   pcdz.Dump(); // call TObject::Dump() to print data members
   pcdz.HoleR=5*GeFiCa::mm; pcdz.HoleH=3*GeFiCa::cm; // add bore hole
   pcdz.HoleTaperW=2*GeFiCa::mm; pcdz.HoleTaperH=2*GeFiCa::mm; // add taper
   TTree *t = pcdz.GetTree(); // create a ROOT tree for quick check
   gStyle->SetPadRightMargin(0.12); // give enough space for color palette
   t->Draw("c2:c1:v","","colz"); // visualize initial setup
}

// show default grid setup of 2D point contact detectors
{
   GeFiCa::PointContactDZ pcdz;
   pcdz.Dump(); // call TObject::Dump() to print data members
   TTree *t = pcdz.GetTree(); // create a ROOT tree for quick check
   t->Draw("c2:c1:v","","colz"); // visualize initial setup
}

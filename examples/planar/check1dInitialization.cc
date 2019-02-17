// show default configuration & grid setup of 1d planar detectors
{
   GeFiCa::Planar1D p1d; // construct a 1D planar detector
   p1d.Dump(); // call TObject::Dump() to print data members

   // create a ROOT tree with the following branches:
   // v: potential
   // e1: electric field in 1st coordintate
   // c1: first coordinate
   // d: depletion flag
   // b: boundary flag
   // i: impurity
   TTree *t = p1d.GetTree();
   t->Scan("v:e1:c1:d:b:i");

   cout<<"Use the following cmd to visualize intial v distribution"<<endl;
   cout<<"t->Draw(\"v:c1\")"<<endl;
}

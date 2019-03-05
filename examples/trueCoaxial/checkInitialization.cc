// show default configuration & grid setup of 1d true coaxial detectors
{
   GeFiCa::TrueCoaxial1D tc1; // construct a 1D true coaxial detector
   tc1.Dump(); // call TObject::Dump() to print data members
   TTree *t = tc1.GetTree(); // create a ROOT tree for quick check
   cout<<"Use the following cmd to visualize intial v distribution"<<endl;
   cout<<"t->Draw(\"v:c1\")"<<endl;
}

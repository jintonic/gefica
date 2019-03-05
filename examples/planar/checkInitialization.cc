// show default configuration & grid setup of 1d planar detectors
{
   GeFiCa::Planar1D p1d; // construct a 1D planar detector
   p1d.Dump(); // call TObject::Dump() to print data members
   TTree *t = p1d.GetTree(); // create a ROOT tree for quick investigation
   cout<<"Use the following cmds for quick investigation"<<endl;
   cout<<"t->Draw(\"v:c1\") or t->Scan(\"v:e:c1:d:b\") etc."<<endl;
}

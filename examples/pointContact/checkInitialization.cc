// show default grid setup of 2D point contact detectors
{
   GeFiCa::PointContactDZ pcdz;
   pcdz.Dump(); // call TObject::Dump() to print data members

   // create a ROOT tree with the following branches:
   // v: potential
   // e1, e2: electric field
   // c1, c2: coordinates
   // d: depletion flag
   // b: boundary flag
   // i: impurity
   TTree *t = pcdz.GetTree();
   t->Draw("c2:c1:v","","colz");
   t->Show(0);
   cout<<"Use the following cmd to visualize intial v distribution"<<endl;
   cout<<"t->Draw(\"v:c1\")"<<endl;
}

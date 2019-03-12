using namespace GeFiCa;
// convert GeFiCa output to fieldgen format for siggen
void g2f(const char *input="ppc.root", const char *output="ppc.dat")
{
   // open input ROOT file
   TFile *fin;
   if (gSystem->Which(".",input)) {
      fin = new TFile(input, "update");
   } else {
      Printf("Cannot open %s", input);
      return;
   }
   // load object from ROOT file
   PointContactDZ *pcdz = (PointContactDZ*) fin->Get("pcdz");
   if (pcdz==0) {
      Printf("Cannot find object in %s", input);
      return;
   }
   TTree *t = pcdz->GetTree();
   double c1, c2, v, e, e1, e2;
   t->SetBranchAddress("potential",&v);
   t->SetBranchAddress("total E  ",&e);
   t->SetBranchAddress("E1            ",&e1);
   t->SetBranchAddress("1st coordinate",&c1);
   t->SetBranchAddress("E2            ",&e2);
   t->SetBranchAddress("2nd coordinate",&c2);
   t->GetEntry(0); double c10=c1; t->GetEntry(1); double c11=c1;
   // dump data to output
   ofstream fout(output);
   fout<<"# xtal_length "<< pcdz->Height/mm<<endl;
   fout<<"# xtal_radius "<<pcdz->Radius/mm<<endl;
   fout<<"# pc_length   "<<pcdz->PointContactH/mm<<endl;
   fout<<"# pc_radius   "<<pcdz->PointContactR/mm<<endl;
   fout<<"# taper_length "<<pcdz->TaperH/mm<<endl;
   fout<<"# wrap_around_radius "<<pcdz->WrapArroundR/mm<<endl;
   fout<<"# ditch_depth "<<pcdz->GrooveH/mm<<endl;
   fout<<"# ditch_thickness "<<pcdz->GrooveW/mm<<endl;
   fout<<"# xtal_grid "<<(c11-c10)/mm<<endl;
   fout<<"# xtal_HV      "<<abs(pcdz->V0-pcdz->V1)/volt<<endl;
   fout<<"# max_iterations "<<pcdz->GetNsor()<<endl;
   fout<<"#"<<endl;
   fout<<"## r (mm), z (mm), V (V),  E (V/cm), E_r (V/cm), E_z (V/cm)"<<endl;
   for (int i=0; i<t->GetEntries(); i++) {
      t->GetEntry(i);
      fout<<c1<<" \t "<<c2<<" \t "<<v<<" \t "<<e<<" \t "<<e1<<" \t "<<e2<<endl;
   }
   fout.close();
   cout<<"data written to "<<output<<endl;
}

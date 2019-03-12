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
   // ROOT to fieldgen format
   pcdz->Export2fieldgen(output);
}

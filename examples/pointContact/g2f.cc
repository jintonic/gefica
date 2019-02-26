using namespace GeFiCa;
/**
 * \file g2f.cc
 * \example g2f.cc
 * \brief Convert GeFiCa output to fieldgen format
 */
void g2f(const char *input="pcdz.root", const char *output="pcdz.dat")
{
   // open input ROOT file
   TFile *fin = 0;
   if (FILE *file = fopen(input,"r")) {
      fclose(file);
      fin = new TFile(input);
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

   // dump object to output
   ofstream fout(output);
   fout.close();
}

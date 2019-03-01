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
   /*
    * fout<<"# height "<< pcdz->Z-pcdz->Z0;        
   fout<<"\n# xtal_radius "<<pcdz->Radius;
   fout<<"\n# pc_length   "<<pcdz->PointContactZ;        
   fout<<"\n# pc_radius   "<<pcdz->PointContactR;         
   fout<<"\n# wrap_around_radius "<<pcdz->WrapArroundR; 
   fout<<"\n# grid size on r "<<pcdz->fdC1p[0];
   fout<<"\n# grid size on z "<<pcdz->fdC2p[0];
   fout<<"\n# impurity_z0  "<<pcdz->fImpurity[0];
   fout<<"\n# xtal_HV      "<<pcdz->V1;
   fout<<"\n# max_iterations "<<pcdz->MaxIterations;
   fout<<"\n# ";
   fout<<"\n## r (mm), z (mm), V (V),  E (V/cm), E_r (V/cm), E_z (V/cm)";
   for (int i=0;i<n;i++) {
      double E=sqrt(fE1[i]*fE1[i]+fE2[i]*fE2[i]);
      fout<<"\n"<<pcdz->fC1[i]<<"  "<<pcdz->fC2[i]<<"  "<<pcdz->fV[i]<<"  "<<pcdz->E<<"  "<<pcdz->fE1[i]<<"  "<<pcdz->fE2[i];
   }
   */
   fout.close();
}

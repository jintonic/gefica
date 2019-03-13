using namespace GeFiCa;
// calculate and save fields of a PPC
void calculateFields(const char *output="ppc.root")
{
   PointContactDZ *ppc = new PointContactDZ(404,501);
   ppc->V1=2.5*kV; ppc->V0=0;

   ppc->Radius=3.45*cm;
   ppc->Height=5.05*cm;

   ppc->PointContactR=1.40*mm;
   ppc->PointContactH=2.10*mm;

   ppc->WrapArroundR=3.45;//1.20*cm;
   ppc->TaperW=0;//0.1*cm;
   ppc->TaperH=0;//0.1*cm;
   ppc->CornerW=0;//0.1*cm;
   ppc->CornerH=0;//0.1*cm;

   ppc->HoleH=0;//2.5*cm;
   ppc->HoleR=0;//0.5*cm;
   ppc->HoleTaperH=0;//0.2*cm;
   ppc->HoleTaperW=0;//0.2*cm;

   ppc->GrooveH=0;//0.1*cm;
   ppc->GrooveW=0;//0.1*cm;

   // x in TF3 -> r, y in TF3 -> z
   TF3 *fid = new TF3("fImpDistr","-0.318e10+0.025e10*y");
   ppc->SetImpurity(fid);

   ppc->Csor=1.995;
   ppc->CalculatePotential(kSOR2);
   
   TFile *file = new TFile(output,"recreate");
   ppc->Write();
   file->Close();
}

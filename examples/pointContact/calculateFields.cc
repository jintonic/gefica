using namespace GeFiCa;
// calculate and save fields of a PPC
void calculateFields(const char *output="ppc.root")
{
   PointContactDZ *ppc = new PointContactDZ(400,501);
   ppc->V0=1.4*kV; ppc->V1=0;

   ppc->Radius=4.00*cm;
   ppc->Height=5.00*cm;

   ppc->PointContactR=0.30*cm;
   ppc->PointContactH=0.10*cm;

   ppc->WrapAroundR=1.20*cm;
   ppc->TaperW=0.1*cm;
   ppc->TaperH=0.1*cm;
   ppc->CornerW=0.1*cm;
   ppc->CornerH=0.1*cm;

   ppc->HoleH=2.5*cm;
   ppc->HoleR=0.5*cm;
   ppc->HoleTaperH=0.2*cm;
   ppc->HoleTaperW=0.2*cm;

   ppc->GrooveH=0.1*cm;
   ppc->GrooveW=0.1*cm;

   // x in TF3 -> r, y in TF3 -> z
   TF3 *fid = new TF3("fImpDistr","0.7e10");
   ppc->SetImpurity(fid);

   ppc->Csor=1.995;
   ppc->CalculatePotential();
   
   TFile *file = new TFile(output,"recreate");
   ppc->Write();
   file->Close();
}

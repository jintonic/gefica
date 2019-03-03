using namespace GeFiCa;
// calculate and save fields of a PPC
void calculateFields(const char *output="ppc.root")
{
   PointContactDZ *fields = new PointContactDZ(692,506);

   fields->Radius=3.45*cm;
   fields->Height=5.05*cm;

   fields->PointContactR=0.14*cm;
   fields->PointContactH=0.21*cm;

   fields->WrapArroundR=2.45*cm;
   fields->TaperW=0.5*cm;
   fields->TaperH=0.5*cm;
   fields->CornerW=0.5*cm;
   fields->CornerH=0.5*cm;

   fields->HoleH=4.0*cm;
   fields->HoleInnerR=0.3*cm;
   fields->HoleOuterR=0.5*cm;

   // x in TF3 -> r, y in TF3 -> z
   TF3 *fid = new TF3("fImpDistr","-0.318e10+0.025e10*y");
   fields->SetImpurity(fid);

   fields->Csor=1.995;
   fields->CalculatePotential(kSOR2);
   
   TFile *file = new TFile(output,"recreate");
   fields->Write();
   file->Close();
}

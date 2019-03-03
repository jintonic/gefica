// calculate and save fields of a PPC
{
   GeFiCa::PointContactDZ *fields = new GeFiCa::PointContactDZ(692,506);
   fields->Radius=3.45*GeFiCa::cm;
   fields->Height=5.05*GeFiCa::cm;
   fields->PointContactR=0.14*GeFiCa::cm;
   fields->PointContactH=0.21*GeFiCa::cm;

   fields->WrapArroundR=2.45*GeFiCa::cm;
   fields->TaperW=0.5*GeFiCa::cm;
   fields->TaperH=0.5*GeFiCa::cm;

   fields->CornerW=0.5*GeFiCa::cm;
   fields->CornerH=0.5*GeFiCa::cm;
   fields->HoleH=4.0*GeFiCa::cm;
   fields->HoleInnerR=0.3*GeFiCa::cm;
   fields->HoleOuterR=0.5*GeFiCa::cm;

   TF3 *im = new TF3("f","0.318e10-0.025e10*y");
   fields->SetImpurity(im);

   fields->CalculatePotential(GeFiCa::kSOR2);
   
   TFile *output = new TFile("ppc.root","recreate");
   fields->Write();
   output->Close();
}

{
   GeFiCa::PointContactDZ *pc=new GeFiCa::PointContactDZ(300,300);
   pc->Radius=3.45*GeFiCa::cm;
   pc->Height=5.05*GeFiCa::cm;
   pc->PointContactR=0.14*GeFiCa::cm;
   pc->PointContactH=0.21*GeFiCa::cm;

   pc->WrapArroundR=2.45*GeFiCa::cm;
   pc->TaperW=0.5*GeFiCa::cm;
   pc->TaperH=0.5*GeFiCa::cm;

   pc->CornerW=0.5*GeFiCa::cm;
   pc->CornerH=0.5*GeFiCa::cm;
   pc->HoleH=4.0*GeFiCa::cm;
   pc->HoleInnerR=0.3*GeFiCa::cm;
   pc->HoleOuterR=0.5*GeFiCa::cm;

   TF3 *im=new TF3("f","-0.318e10+0.025e10*y");
   pc->SetImpurity(im);

   pc->CalculatePotential(GeFiCa::kSOR2);
   
   TFile * outfile=new TFile("rcpc.root","UPDATE");
   pc->Write("rcpc");
   outfile->Write();

}

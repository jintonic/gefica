void compare2fieldgen()
{
   double r,z,p;

   TH2F *hmjd = new TH2F("hmjd","",346,0,346,506,0,506);
   ifstream ff("ev.dat");
   for (int i = 0; i < 346; i++) {
     for (int j = 0; j < 506; j++) {
        ff>>r>>z>>p;
        hmjd->SetBinContent(i,j,p);
     } 
   }

   TH2F *hgfc = new TH2F("hgfc","",346,0,346,506,0,506);
   ifstream fg("halfP.txt");
   for (int i = 0; i < 346; i++) {
     for (int j = 0; j < 506; j++) {
        fg>>r>>z>>p;
        hgfc->SetBinContent(i,j,p);
     } 
   }

   hmjd->Draw("colz");
   TCanvas *can = new TCanvas;
   hgfc->Draw("colz");
}

void generateField()
{
   GeFiCa::halfPointContactRZ *ppc = new GeFiCa::halfPointContactRZ(1036,506);
   ppc->Radius=3.45;
   ppc->ZUpperBound=5.05;
   ppc->PointR=0.135;
   ppc->PointDepth=1.05;

   ppc->MaxIterations=1e6;
   ppc->Precision=1e-8;
   ppc->Csor=1.996;
   ppc->V0=0*GeFiCa::volt;
   ppc->V1=-2500*GeFiCa::volt;
   ppc->Impurity="-0.318e10+0*y";

   ppc->CalculateField(GeFiCa::kSOR2);
   ppc->SaveField("half.root");
}

void r2t()
{
   GeFiCa::halfPointContactRZ *ppc = new GeFiCa::halfPointContactRZ(1036,506);
   ppc->LoadField("half.root");

   ofstream fo("halfP.txt");
   for (int i = 0; i < 346; i++) {
     for (int j = 0; j < 506; j++) {
        fo<<i*ppc->Radius/346<<"\t"<<j*ppc->ZUpperBound/506<<"\t"
           <<ppc->GetPotential(i*ppc->Radius/346,j*ppc->ZUpperBound/506)<<endl;
     } 
   }
   fo.close();
}

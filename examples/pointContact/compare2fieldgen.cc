using namespace GeFiCa;
void compare2fieldgen(const char *fieldgenOutput="fieldgen/fields/p1/ev.dat")
{
   PointContact detector; // create a point contact detector
   // same as fieldgen/config_files/p1_new.config
   detector.Bias[0]=0*kV; // bias on point contact
   detector.Bias[1]=2.5*kV; // bias on outer contact
   detector.Radius=34.5*mm; detector.Height=50.5*mm;
   detector.PointContactR=1.4*mm; detector.PointContactH=2.1*mm;
   // what's set in fieldgen is space charge instead of impurity density
   detector.BottomImpurity=0.318e10/cm3; detector.TopImpurity=0.19175e10/cm3;
   // no fancy feature
   detector.TaperW=0; detector.TaperH=0;
   detector.CornerW=0; detector.CornerH=0;
   detector.BoreH=0; detector.BoreR=0;
   detector.BoreTaperH=0; detector.BoreTaperW=0;
   detector.GrooveH=0; detector.GrooveW=0;

   RhoZ grid(690,506); // grid step length is about 0.1 mm
   grid.GetBoundaryConditionFrom(detector);
   grid.RelaxationFactor=1.994; // fastest
   grid.SuccessiveOverRelax();

   // get fieldgen data
   ifstream fieldgen(fieldgenOutput); if (!fieldgen.is_open()) exit(-1);
   string str; for (int i=0; i<36; i++) getline(fieldgen,str); // skip header

   // compare
   ofstream output("gVSf.txt");
   double r,z,v,e,er,ez,V,E,Er,Ez;
   while (fieldgen>>r>>z>>v>>e>>er>>ez) {
      V=grid.GetV(r*mm,z*mm);
      Er=grid.GetE1(r*mm,z*mm);
      Ez=grid.GetE2(r*mm,z*mm);
      E=sqrt(Er*Er+Ez*Ez);
      output<<r<<"  "<<z<<"  "<<V<<"  "<<v-V<<"  "<<Er<<"  "
         <<Ez<<"  "<<E<<"  "<<er-Er<<"  "<<ez-Ez<<"  "<<e-E<<endl;
   }
   fieldgen.close();
   output.close();

   gROOT->SetStyle("GeFiCa");
   gStyle->SetTitleOffset(0.8,"Y");
   gStyle->SetPadRightMargin(0.15);
   gStyle->SetPadLeftMargin(0.08);

   // draw fieldgen result
   TCanvas *c1 = new TCanvas("c1", "fieldgen");
   TTree *tf = new TTree("tf","tf");
   tf->ReadFile(fieldgenOutput, "r:z:v:e:er:ez");
   tf->Draw("z:r:v","","colz");

   // draw GeFiCa result
   TCanvas *c2 = new TCanvas("c2", "gefica");
   TTree *tg = new TTree("tg","tg");
   tg->ReadFile("gVSf.txt", "r:z:v:d:e1:e2:e:de1:de2:de");
   tg->Draw("z:r:v","","colz");

   // draw difference
   TCanvas *c3 = new TCanvas("c3", "fieldgen - gefica");
   tg->Draw("z:r:d","","colz"); // potential difference

   TCanvas *c4 = new TCanvas("c4","(fieldgen - gefica)/fieldgen*100",600,600);
   c4->SetLeftMargin(0.1); c4->SetRightMargin(0.2);
   tg->Draw("z:r:de/e*100","(r>1.4 || z>2.3)","goff"); // field difference
   TGraph2D *gg = new TGraph2D(tg->GetEntries(),
         tg->GetV2(), tg->GetV1(), tg->GetV3());
   gg->SetName("gg"); gg->SetNpx(500); gg->SetNpy(500); // fine bin histogram
   TH2D *he = gg->GetHistogram();
   he->SetTitle(";Radius [cm];Height [cm];"
         "100 #times (E_{fieldgen}-E_{GeFiCa})/E_{fieldgen}");
   he->GetZaxis()->CenterTitle();
   he->GetZaxis()->SetTitleOffset(1.5);
   he->GetYaxis()->SetTitleOffset(0.9);
   he->Draw("colz2");

   TBox *box = new TBox(0.05,0,1.42,2.3);
   box->SetFillColor(kWhite);
   box->Draw();

   c4->Print("gVSf.png");
}

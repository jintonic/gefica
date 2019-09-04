using namespace GeFiCa;
//______________________________________________________________________________
// config: pancake PPC with large PC
void compare2pancake()
{
   PointContact detector;
   double nRho=-0.318e10/cm3; // space charge density of a p-type detector
   detector.SetAverageImpurity(-nRho);
   detector.Radius=21*mm; detector.Height=5*mm; // pancake
   detector.PointContactR=20*mm; detector.PointContactH=0*mm; // large PC
   detector.Bias[0]=0*kV; // bias on point contact
   detector.Bias[1]=100*volt; // bias on outer contact
   // no fancy feature
   detector.TaperW=0; detector.TaperH=0;
   detector.CornerW=0; detector.CornerH=0;
   detector.BoreH=0; detector.BoreR=0;
   detector.BoreTaperH=0; detector.BoreTaperW=0;
   detector.GrooveH=0; detector.GrooveW=0;

   const int nr=420, nz = 51; // number of grid points
   RhoZ grid(nr,nz); // grid step length is about 0.1 mm
   grid.GetBoundaryConditionFrom(detector);
   grid.RelaxationFactor=1.994; // fastest
   grid.SuccessiveOverRelax();

   gROOT->SetStyle("GeFiCa");
   gStyle->SetLabelFont(133,"XYZ");
   gStyle->SetTitleFont(133,"XYZ");
   gStyle->SetLabelSize(18,"XYZ");
   gStyle->SetTitleSize(18,"XYZ");
   gStyle->SetPadRightMargin(0.005);
   gStyle->SetPadLeftMargin(0.09);
   gStyle->SetPadBottomMargin(0.15);

   // analytic solution
   TCanvas *can1 = new TCanvas("can1","can1",600,300);
   double z[nz], va[nz], vg[nz];
   double a =-nRho*Qe/epsilon;
   double b = detector.Bias[1]/detector.Height - a/2*detector.Height;
   for (int i=0; i<nz; i++) {
      z[i]=i/10.;
      va[i]=a*z[i]*mm*z[i]*mm/2 + b*z[i]*mm;
      vg[i]=grid.GetV(0, z[i]*mm);
   }
   TGraph *ga = new TGraph(nz,z,va);
   TGraph *gg = new TGraph(nz,z,vg);
   ga->SetTitle(";Axial position at radius=0 [mm];Voltage [volt]");
   ga->GetYaxis()->SetTitleOffset(0.7);
   ga->SetLineColor(kRed);
   gg->SetMarkerColor(kBlue);
   ga->Draw("al");
   gg->Draw("p");

   TLegend *l = new TLegend(0.3,0.7,0.5,0.85);
   l->AddEntry(ga,"Analytic","l");
   l->AddEntry(gg,"GeFiCa","p");
   l->Draw();

   can1->Print("pancakeVZ.png");

   // GeFiCa result
   gStyle->SetPadTopMargin(0.08);
   gStyle->SetPadBottomMargin(0.32);
   gStyle->SetPadLeftMargin(0.08);
   gStyle->SetPadRightMargin(0.1);
   gStyle->SetTitleOffset(-0.13,"Z");
   TCanvas *can2 = new TCanvas("can2","can2",600,140);
   TTree *t = grid.GetTree();
   t->Draw("c2*10:c1*10:v","","goff");
   TGraph2D *g2 = new TGraph2D(t->GetEntries(),
         t->GetV2(), t->GetV1(), t->GetV3());
   g2->SetName(""); g2->SetNpx(500); g2->SetNpy(500); // fine bin histogram
   TH2D *he = g2->GetHistogram();
   he->SetTitle(";Radius [mm];Height [mm];V [volt]");
   he->GetYaxis()->SetTitleOffset(0.23);
   he->GetYaxis()->SetNdivisions(5);
   he->GetZaxis()->CenterTitle();
   he->GetZaxis()->SetNdivisions(5);
   he->Draw("colz");
   can2->Print("pancake.png");
}
//______________________________________________________________________________
// config: long & thin PPC with tall PC
void compare2match()
{
   PointContact detector;
   double nRho=-0.318e10/cm3; // space charge density of a p-type detector
   detector.SetAverageImpurity(-nRho);
   detector.Radius=5*mm; detector.Height=50*mm; // match shape
   detector.PointContactR=2*mm; detector.PointContactH=49*mm; // tall PC
   detector.Bias[0]=0*kV; // bias on point contact
   detector.Bias[1]=100*volt; // bias on outer contact
   // no fancy feature
   detector.TaperW=0; detector.TaperH=0;
   detector.CornerW=0; detector.CornerH=0;
   detector.BoreH=0; detector.BoreR=0;
   detector.BoreTaperH=0; detector.BoreTaperW=0;
   detector.GrooveH=0; detector.GrooveW=0;

   const int nr=100, nz = 51; // number of grid points
   RhoZ grid(nr,nz); // grid step length is about 0.1 mm
   grid.GetBoundaryConditionFrom(detector);
   grid.RelaxationFactor=1.994; // fastest
   grid.SuccessiveOverRelax();

   gROOT->SetStyle("GeFiCa");
   gStyle->SetLabelFont(133,"XYZ");
   gStyle->SetTitleFont(133,"XYZ");
   gStyle->SetLabelSize(18,"XYZ");
   gStyle->SetTitleSize(18,"XYZ");
   gStyle->SetPadRightMargin(0.005);
   gStyle->SetPadLeftMargin(0.22);
   gStyle->SetPadTopMargin(0.018);
   gStyle->SetPadBottomMargin(0.07);

   // analytic solution
   TCanvas *can3 = new TCanvas("can3","can3",230,600);
   double r[nr], va[nr], vg[nr];
   double a =-nRho*Qe/epsilon;
   double b1= (detector.Bias[1] - a*(detector.Radius*detector.Radius
            -detector.PointContactR*detector.PointContactR)/4)
      /log(detector.Radius/detector.PointContactR);
   double b2= detector.Bias[1]*log(detector.PointContactR)
      - a*(detector.Radius*detector.Radius*log(detector.PointContactR)
            -detector.PointContactR*detector.PointContactR
            *log(detector.Radius))/4;
   b2/=log(detector.PointContactR)-log(detector.Radius);
   for (int i=0; i<nr/2-20; i++) {
      r[i]=i/10.+2;
      va[i]=a*r[i]*mm*r[i]*mm/4 + b1*log(r[i]*mm) + b2;
      vg[i]=grid.GetV(r[i]*mm,5*mm);
   }
   TGraph *ga = new TGraph(nr/2-20,r,va);
   TGraph *gg = new TGraph(nr/2-20,r,vg);
   ga->SetTitle(";R position@z=5mm [mm];Voltage [volt]");
   ga->GetYaxis()->SetTitleOffset(2.9);
   ga->SetLineColor(kRed);
   gg->SetMarkerColor(kBlue);
   ga->Draw("al");
   gg->Draw("p");

   TLegend *l = new TLegend(0.5,0.2,0.99,0.35);
   l->SetTextFont(133);
   l->SetTextSize(18);
   l->AddEntry(ga,"Analytic","l");
   l->AddEntry(gg,"GeFiCa","p");
   l->Draw();

   can3->Print("matchVR.png");

   // GeFiCa result
   gStyle->SetPadTopMargin(0.02);
   gStyle->SetPadBottomMargin(0.07);
   gStyle->SetPadLeftMargin(0.19);
   gStyle->SetPadRightMargin(0.3);
   gStyle->SetTitleOffset(-1.5,"Z");
   TCanvas *can4 = new TCanvas("can4","can4",200,600);
   TTree *t = grid.GetTree();
   t->Draw("c2*10:c1*10:v","","goff");
   TGraph2D *g2 = new TGraph2D(t->GetEntries(),
         t->GetV2(), t->GetV1(), t->GetV3());
   g2->SetName(""); g2->SetNpx(500); g2->SetNpy(500); // fine bin histogram
   TH2D *he = g2->GetHistogram();
   he->SetTitle(";Radius [mm];Height [mm];V [volt]");
   he->GetYaxis()->SetTitleOffset(2.7);
   he->GetZaxis()->CenterTitle();
   he->GetXaxis()->SetNdivisions(3);
   he->Draw("colz");

   gPad->Modified(); gPad->Update();
   TPaletteAxis *palette = (TPaletteAxis*)
      he->GetListOfFunctions()->FindObject("palette");
   palette->SetX1NDC(0.71);
   palette->SetX2NDC(0.86);

   can4->Print("match.png");
}
//______________________________________________________________________________
//
void compare2analytic()
{
   compare2pancake();
   compare2match();
}

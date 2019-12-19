// definition of necessary units
static const double cm=1;
static const double cm3=cm*cm*cm;
static const double volt=1;
static const double C=1; // Coulomb
static const double Qe=-1.6e-19*C; // eletron charge
static const double epsilon0=8.854187817e-14*C/volt/cm; // vacuum permittivity
// https://link.springer.com/chapter/10.1007/10832182_519
static const double epsilonGe=15.8; // Ge dielectric constant
//______________________________________________________________________________
// V"(x)=a, https://www.wolframalpha.com/input/?i=V%27%27(x)%3Da
double V(double *coordinates, double *parameters)
{
	double x = coordinates[0];  // there is no y and z dependence
	double x0= 0*cm;            // lower electrode
	double x1= parameters[0];   // upper electrode
	double v0= 0*volt;          // lower voltage
	double v1= parameters[1];   // upper voltage
	double rho=parameters[2]*Qe;// space charge density [C/cm3]

	double a =-rho/epsilon0/epsilonGe;
	double c2= (v1-v0)/(x1-x0) - a/2*(x1+x0);
	double c1= (v0*x1-v1*x0)/(x1-x0) + a/2*x0*x1;
	return a*x*x/2 + c2*x + c1;
}
//______________________________________________________________________________
// search for depletion voltage of a planar detector given impurity & thickness
double GetVdep(double impurity, double thickness, double vupper)
{
	if (impurity==0) return 0; // nothing to deplete

	TF1 *det=new TF1("det", V, 0, thickness, 3); // potential distr.
	double bias, vlower=0*volt; // range of search
	while (abs(vupper-vlower)>1e-3*volt) { // binary search
		bias=(vupper+vlower)/2; // bias voltage
		det->SetParameters(thickness, bias, impurity);
		if (det->Derivative(0)*det->Derivative(thickness)<0) vlower=bias;
		else vupper=bias;
	}
	return bias;
}
//______________________________________________________________________________
// draw results
void getVd(int type=1, double thickness=1*cm)
{
	const int n=10; // number of points
	double impurity[n]={1e9/cm3, 2e9/cm3, 4e9/cm3, 8e9/cm3, 1e10/cm3,
		2e10/cm3, 4e10/cm3, 8e10/cm3, 1e11/cm3, 1.2e11/cm3}; // p-type
	if (type==-1) for (int i=0; i<n; i++) impurity[i]*=type; // n-type
	double absi[n], vdep[n];
	for (int i=0; i<n; i++) {
		vdep[i] = GetVdep(impurity[i], thickness, type*-2e4*volt);
		absi[i] = abs(impurity[i]);
	}

	gROOT->SetStyle("Plain"); // pick up a good default drawing style
	// modify the default style
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(132);
	gStyle->SetLabelFont(132,"XY");
	gStyle->SetTitleFont(132,"XY");
	gStyle->SetLabelSize(0.05,"XY");
	gStyle->SetTitleSize(0.05,"XY");
	gStyle->SetTitleOffset(1.1,"XY");
	gStyle->SetPadRightMargin(0.01);
	gStyle->SetPadLeftMargin(0.11);
	gStyle->SetPadTopMargin(0.01);
	gStyle->SetPadBottomMargin(0.12);

	// Vdep VS impurity
	TGraph *g = new TGraph(n,absi,vdep);
	g->SetTitle(";Net Impurity [cm^{-3}];Depletion Voltage [V]");
	g->Draw("apc");
	TText *t1 = new TText(2e9, 6000, Form(
				"%.0f cm thick planar detector", thickness/cm));
	t1->Draw();
	gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
	gPad->Print("depleted.png");

	// voltage VS thickness
	TCanvas *c = new TCanvas;
	double bias = -1000*volt; // p-type
	if(type==-1) bias = 1000*volt; // n-type

	// over depleted
	TF1 *fvo=new TF1("fvo", V, 0, thickness, 3);
	fvo->SetParameters(thickness, vdep[6]+bias, impurity[6]);
	fvo->SetTitle(";Thickness [cm]; Voltage [V]");
	fvo->SetLineColor(kMagenta);
	fvo->SetLineStyle(2);
	fvo->GetXaxis()->SetTitleOffset(1.1);
	fvo->Draw();

	// depleted
	TF1 *fvd=new TF1("fvd", V, 0, thickness, 3);
	fvd->SetParameters(thickness, vdep[6], impurity[6]);
	fvd->Draw("same");

	// undepleted
	TF1 *fvu=new TF1("fvu", V, 0, thickness, 3);
	fvu->SetParameters(thickness, bias, impurity[6]);
	fvu->SetLineColor(kRed);
	fvu->SetLineStyle(3);
	fvu->Draw("same");

	// potential due to bias alone
	TF1 *fvb=new TF1("fvb", V, 0, thickness, 3);
	fvb->SetParameters(thickness, bias, 0/cm3);
	fvb->SetLineColor(kBlue);
	fvb->SetLineStyle(4);
	fvb->Draw("same");

	// potential due to space charges alone
	TF1 *fvc=new TF1("fvc", V, 0, thickness, 3);
	fvc->SetParameters(thickness, 0*volt, impurity[6]);
	fvc->SetLineColor(kGreen);
	fvc->SetLineStyle(5);
	fvc->Draw("same");

	if (type==-1) {
		// draw lines and texts for n-type
		TLine *l1 = new TLine(0,bias/volt,thickness/cm,bias/volt);
		l1->SetLineStyle(kDashed); l1->Draw();
		TLine *l2 = new TLine(0,vdep[6]/volt,thickness/cm,vdep[6]/volt);
		l2->SetLineStyle(kDashed); l2->Draw();
		TLine *l3=new TLine(0,(vdep[6]+bias)/volt,thickness/cm,(vdep[6]+bias)/volt);
		l3->SetLineStyle(kDashed); l3->Draw();
		TText *t2 = new TText(0.04,bias/volt+20,Form("%.0f V", bias/volt));
		t2->SetTextFont(132); 
		t2->Draw();
		TText *t3 = new TText(0.04,vdep[6]/volt+20, Form("%.0f V", vdep[6]/volt));
		t3->SetTextFont(132); 
		t3->Draw();
		TText *t4 = new TText(0.04,(vdep[6]+bias)/volt+20,
				Form("%.0f V", (vdep[6]+bias)/volt));
		t4->Draw();
		TText *t5 = new TText(0.5, 2800, "over depleted");
		t5->SetTextFont(132); 
		t5->SetTextColor(kMagenta); t5->Draw();
		TText *t6 = new TText(0.8, 2310, "just depleted");
		t6->SetTextFont(132); 
		t6->Draw();
		TText *t7 = new TText(0.65, 1250, "undepleted");
		t7->SetTextFont(132); 
		t7->SetTextColor(kRed); t7->Draw();
		TText *t8 = new TText(0.65, 800, "bias alone");
		t8->SetTextFont(132); 
		t8->SetTextColor(kBlue); t8->Draw();
		TText *t9 = new TText(0.55, 200, "space charges alone");
		t9->SetTextFont(132); 
		t9->SetTextColor(kGreen); t9->Draw();
		TLatex *t10 = new TLatex(0.1, 2700,
				Form("Impurity: %.0e/cm^{3}",impurity[6]/cm3));
		t10->SetTextFont(132); 
		t10->Draw();
	} else {
		// draw lines and texts for p-type
		TLine *l0 = new TLine(0,0,thickness/cm,0);
		l0->SetLineStyle(kDashed); l0->Draw();
		TLine *l1 = new TLine(0,bias/volt,thickness/cm,bias/volt);
		l1->SetLineStyle(kDashed); l1->Draw();
		TLine *l2 = new TLine(0,vdep[6]/volt,thickness/cm,vdep[6]/volt);
		l2->SetLineStyle(kDashed); l2->Draw();
		TLine *l3=new TLine(0,(vdep[6]+bias)/volt,thickness/cm,(vdep[6]+bias)/volt);
		l3->SetLineStyle(kDashed); l3->Draw();
		TText *t2 = new TText(0.04,bias/volt+20,Form("%.0f V", bias/volt));
		t2->SetTextFont(132); 
		t2->Draw();
		TText *t3 = new TText(0.04,vdep[6]/volt+20, Form("%.0f V", vdep[6]/volt));
		t3->SetTextFont(132); 
		t3->Draw();
		TText *t4 = new TText(0.04,(vdep[6]+bias)/volt+20,
				Form("%.0f V", (vdep[6]+bias)/volt));
		t4->SetTextFont(132); 
		t4->Draw();
		TText *t5 = new TText(0.5, -3100, "over depleted");
		t5->SetTextFont(132); 
		t5->SetTextColor(kMagenta); t5->Draw();
		TText *t6 = new TText(0.76, -2500, "just depleted");
		t6->SetTextFont(132); 
		t6->Draw();
		TText *t7 = new TText(0.65, -1400, "undepleted");
		t7->SetTextFont(132); 
		t7->SetTextColor(kRed); t7->Draw();
		TText *t8 = new TText(0.6, -900, "bias alone");
		t8->SetTextFont(132); 
		t8->SetTextColor(kBlue); t8->Draw();
		TText *t9 = new TText(0.5, -200, "space charges alone");
		t9->SetTextFont(132); 
		t9->SetTextColor(kGreen); t9->Draw();
		TLatex *t10 = new TLatex(0.1, -2700,
				Form("Impurity: %.0e/cm^{3}",impurity[6]/cm3));
		t10->SetTextFont(132); 
		t10->Draw();
	}
   c->Print("undepleted.png");
}

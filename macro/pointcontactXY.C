using namespace GeFiCa;

void pointcontactXY()
{
	PointContactXY *detector = new PointContactXY(101,101);
	detector->MaxIterations=1e5;
	detector->Csor=1;
//	detector->SetImpurity(1e10/cm3);
	detector->CalculateField(EMethod::kSOR2);
	detector->SaveField("pcXY.root");
	
	TCanvas *c = new TCanvas;
	TChain *t1 = new TChain("t");
	t1->Add("pcXY.root");
	t1->Draw("p:c1:c2");


}
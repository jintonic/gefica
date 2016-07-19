using namespace GEFICA;
void test()
{
   idntnshape * field=new idntnshape(101,10,10);
   field->Csor=1.1;
   double stepLength=0.1;
   field->Create(stepLength);
   field->SetVoltage(2000,1); // field->SetVoltage(2000);
   field->MaxIterations=100000; // field->MaxIterations=1000;
   TF1 *f = new TF1("f","pol1", 0, 10);

   f->SetParameters(1e10, 0);

   //field->SetImpurity(f);
   bool converged=field->Iterate();
   field->Save("out.root");
   if (!converged) {
      field->Load("out.root");
      field->Iterate();
      field->Save("out.root");
      cout<<"1"<<endl;
   }
}

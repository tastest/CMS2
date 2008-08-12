{

  gROOT->SetBatch(kTRUE);

  gROOT->LoadMacro("../Tools/histtools.C");

  TCanvas *canvas = new TCanvas("stacks","stacks",1200.,800.);
  TCanvas *canvas2 = new TCanvas("stack","stack",1200.,800.);

  canvas->Divide(6,4);

  const unsigned int allBuckets = 20;
  char *suffix[allBuckets+1];
  suffix[0]  = "mpmpmp";
  suffix[1]  = "mpmpmm";
  suffix[2]  = "mpmmmm";
  suffix[3]  = "mmmmmm";
  suffix[4]  = "mpmpep";
  suffix[5]  = "mpmpem";
  suffix[6]  = "mpmmep";
  suffix[7]  = "mpmmem";
  suffix[8]  = "mmmmep";
  suffix[9]  = "mmmmem";
  suffix[10] = "mpepep";
  suffix[11] = "mpepem";
  suffix[12] = "mpemem";
  suffix[13] = "mmepep";
  suffix[14] = "mmepem";
  suffix[15] = "mmemem";
  suffix[16] = "epepep";
  suffix[17] = "epepem";
  suffix[18] = "epemem";
  suffix[19] = "ememem";
  suffix[20] = "all";

  const char* histname = "hNjetsBothLeptonsVeto";

  for ( unsigned int bucket = 0;
	bucket < allBuckets+1;
	++bucket ) {
    canvas->cd(bucket+1);
    THStack *stack = hist::stack(Form("st_%s_%s",histname,suffix[bucket]),Form("._%s_%s",histname,suffix[bucket]));
    stack->Draw("HIST");
    TLegend *leg = hist::legend(stack,"lpf",0,0);
    leg->SetBorderSize(1);
    leg->SetFillColor(0);
    leg->Draw("SAME");
    canvas->Update();

    canvas2->Clear();
    canvas2->cd();
    stack->Draw("HIST");
    leg->Draw("SAME");
    canvas2->Update();
    canvas2->Print(Form("stack_%s_%s.png",histname,suffix[bucket]));

  }

  canvas->Print(Form("stacks_%s.png",histname));

  gROOT->SetBatch(kFALSE);

}

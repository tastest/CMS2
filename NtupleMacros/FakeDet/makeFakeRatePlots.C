TCanvas canvas("canvas","canvas",900,900);

void makeFakeRatePlots() {
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);

  TFile *output = new TFile("fakeRates.root","RECREATE");

  vector<TString> samples;

  samples.push_back(TString("fakeRatesFull"));

  for ( unsigned int sample = 0;
	sample < samples.size();
	++sample ) {
    TString filename(samples[sample]);
    filename.Append(".root");
    TFile *file = new TFile(filename);
    if ( file->IsZombie() ) {
      cout << "File: " << filename << " cannot be opened and is skipped!" << endl;
    } else {
      makeFractionPlots(file,		      
			output,
			samples[sample],
			"num_ell",
			"den_ele",
			"#eta",    
			"p_{T} [GeV]",    
			"Fake Rate");

      makeFractionPlots(file,		      
			output,
			samples[sample],
			"num_wo_leading_ell",
			"den_wo_leading_ele",
			"#eta",    
			"p_{T} [GeV]",    
			"Fake Rate");

      makeFractionPlots(file,		      
			output,
			samples[sample],
			"num_wo_second_leading_ell",
			"den_wo_second_leading_ele",
			"#eta",    
			"p_{T} [GeV]",    
			"Fake Rate");
      makeFractionPlots(file,		      
			output,
			samples[sample],
			"num_elt",
			"den_ele",
			"#eta",    
			"p_{T} [GeV]",    
			"Fake Rate");

      makeFractionPlots(file,		      
			output,
			samples[sample],
			"num_wo_leading_elt",
			"den_wo_leading_ele",
			"#eta",    
			"p_{T} [GeV]",    
			"Fake Rate");

      makeFractionPlots(file,		      
			output,
			samples[sample],
			"num_wo_second_leading_elt",
			"den_wo_second_leading_ele",
			"#eta",    
			"p_{T} [GeV]",    
			"Fake Rate");
    }
  }

  output->Write();

}

void makeFractionPlots(TFile *file,
		       TFile *output,
		       TString &sample,
		       TString &numerator,
		       TString &denominator,
		       TString &x_axis_title,    
		       TString &y_axis_title,    
		       TString &z_axis_title) {

  canvas.Clear();
  canvas.Divide(2,2);

  TH2F *num = (TH2F*)file->Get(numerator);
  TH2F *den = (TH2F*)file->Get(denominator);
  TH2F *fake_rate = (TH2F*)num->Clone();
  fake_rate->Reset();
  fake_rate->Divide(num,den,1.,1.,"B");
  TString name = "fakeRateTemplate_";
  name.Append(numerator);
  name.Append("_");
  name.Append(sample);
  name.ReplaceAll("num_","");
  TString title = "Fake Rate ";
  title.Append(numerator);
  title.Append(" ");
  title.Append(sample);
  title.ReplaceAll("num_","");
  fake_rate->SetName(name);
  fake_rate->SetTitle(title);
  fake_rate->SetXTitle(x_axis_title);
  fake_rate->SetYTitle(y_axis_title);
  fake_rate->SetZTitle(z_axis_title);
  fake_rate->SetMarkerStyle(20);
  fake_rate->SetMarkerColor(kBlack);
  fake_rate->SetDirectory(output);
  fake_rate->SetMinimum(0.0);
  fake_rate->SetMaximum(1.0);
  fake_rate->GetXaxis()->SetLabelSize(0.03);
  fake_rate->GetXaxis()->SetTitleSize(0.03);
  fake_rate->GetXaxis()->SetTitleOffset(2.);
  fake_rate->GetYaxis()->SetLabelSize(0.03);
  fake_rate->GetYaxis()->SetTitleSize(0.03);
  fake_rate->GetYaxis()->SetTitleOffset(3.);
  fake_rate->GetZaxis()->SetLabelSize(0.04);
  fake_rate->GetZaxis()->SetTitleSize(0.04);
  fake_rate->GetZaxis()->SetTitleOffset(1.35);

  TH2F *fake_rate_error = (TH2F*)fake_rate->Clone();
  fake_rate_error->Reset();
  TString name = "fakeRateTemplateError_";
  name.Append(numerator);
  name.Append("_");
  name.Append(sample);
  name.ReplaceAll("num_","");
  TString title = "Fake Rate error ";
  title.Append(numerator);
  title.Append(" ");
  title.Append(sample);
  title.ReplaceAll("num_","");
  fake_rate_error->SetName(name);
  fake_rate_error->SetTitle(title);
  fake_rate_error->SetXTitle(x_axis_title);
  fake_rate_error->SetYTitle(y_axis_title);
  fake_rate_error->SetZTitle(z_axis_title);
  fake_rate_error->SetMarkerStyle(20);
  fake_rate_error->SetMarkerColor(kBlack);
  fake_rate_error->SetDirectory(output);
  fake_rate_error->SetMinimum(0.0);
  fake_rate_error->SetMaximum(1.0);
  for ( unsigned int x = 0;
	x <= fake_rate->GetNbinsX();
	++x) {
    for ( unsigned int y = 0;
	  y <= fake_rate->GetNbinsY();
	  ++y) {
      fake_rate_error->SetBinContent(x,y,fake_rate->GetBinError(x,y));
    }
  }

  TH1D *fake_rate_project_x_num = num->ProjectionX("1",-1,-1,"e");
  TH1D *fake_rate_project_x_den = den->ProjectionX("2",-1,-1,"e");
  TH1D *fake_rate_project_x = dynamic_cast<TH1D*>fake_rate_project_x_num->Clone("clone_1");
  fake_rate_project_x->Reset();
  fake_rate_project_x->Divide(fake_rate_project_x_num,fake_rate_project_x_den,1.,1.,"B");
  name = "fakeRateTemplate_";
  name.Append(numerator);
  name.Append("_");
  name.Append(sample);
  name.Append("_eta");
  name.ReplaceAll("num_","");
  title = "Fake Rate ";
  title.Append(numerator);
  title.Append(" ");
  title.Append(sample);
  title.Append(" eta");
  title.ReplaceAll("num_","");
  fake_rate_project_x->SetName(name);
  fake_rate_project_x->SetTitle(title);
  fake_rate_project_x->SetXTitle(x_axis_title);
  fake_rate_project_x->SetYTitle(z_axis_title);
  fake_rate_project_x->SetMarkerStyle(20);
  fake_rate_project_x->SetMarkerColor(kRed);
  fake_rate_project_x->SetDirectory(output);
  fake_rate_project_x->SetMinimum(0.0);
  fake_rate_project_x->SetMaximum(0.5);

  TH1D *fake_rate_project_y_num = num->ProjectionY("3",-1,-1,"e");
  TH1D *fake_rate_project_y_den = den->ProjectionY("4",-1,-1,"e");
  TH1D *fake_rate_project_y = dynamic_cast<TH1D*>fake_rate_project_y_num->Clone("clone_2");
  fake_rate_project_y->Reset();
  fake_rate_project_y->Divide(fake_rate_project_y_num,fake_rate_project_y_den,1.,1.,"B");
  name = "fakeRateTemplate_";
  name.Append(numerator);
  name.Append("_");
  name.Append(sample);
  name.Append("_pt");
  name.ReplaceAll("num_","");
  title = "Fake Rate ";
  title.Append(numerator);
  title.Append(" ");
  title.Append(sample);
  title.Append(" pt");
  title.ReplaceAll("num_","");
  fake_rate_project_y->SetName(name);
  fake_rate_project_y->SetTitle(title);
  fake_rate_project_y->SetXTitle(y_axis_title);
  fake_rate_project_y->SetYTitle(z_axis_title);
  fake_rate_project_y->SetMarkerStyle(20);
  fake_rate_project_y->SetMarkerColor(kRed);
  fake_rate_project_y->SetDirectory(output);
  fake_rate_project_y->SetMinimum(0.0);
  fake_rate_project_y->SetMaximum(0.5);

  TPad *pad1 = canvas.cd(1);
  fake_rate->Draw("LEGO2");
  TPaveStats *fake_rate_stats = (TPaveStats*)(fake_rate->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats != 0 ) {
    fake_rate_stats->SetX1NDC(0.7);
    fake_rate_stats->SetY1NDC(0.75);
    fake_rate_stats->SetX2NDC(0.9);
    fake_rate_stats->SetY2NDC(0.99);
  }
  pad1->Update();

  TPad *pad2 = canvas.cd(2);
  fake_rate_error->Draw("LEGO2");
  TPaveStats *fake_rate_error_stats = (TPaveStats*)(fake_rate_error->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_error_stats != 0 ) {
    fake_rate_error_stats->SetX1NDC(0.7);
    fake_rate_error_stats->SetY1NDC(0.75);
    fake_rate_error_stats->SetX2NDC(0.9);
    fake_rate_error_stats->SetY2NDC(0.99);
  }
  pad2->Update();

  TPad *pad3 = canvas.cd(3);
  fake_rate_project_x->Draw("PE");
  TPaveStats *fake_rate_project_x_stats = (TPaveStats*)(fake_rate_project_x->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_project_x_stats != 0 ) {
    fake_rate_project_x_stats->SetX1NDC(0.7);
    fake_rate_project_x_stats->SetY1NDC(0.75);
    fake_rate_project_x_stats->SetX2NDC(0.9);
    fake_rate_project_x_stats->SetY2NDC(0.99);
  }
  pad3->SetLeftMargin(0.15);
  pad3->SetRightMargin(0.1);
  pad3->SetBottomMargin(0.17);
  pad3->SetTopMargin(0.1);
  pad3->Update();

  TPad *pad4 = canvas.cd(4);
  fake_rate_project_y->Draw("PE");
  TPaveStats *fake_rate_project_y_stats = (TPaveStats*)(fake_rate_project_y->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_project_y_stats != 0 ) {
    fake_rate_project_y_stats->SetX1NDC(0.7);
    fake_rate_project_y_stats->SetY1NDC(0.75);
    fake_rate_project_y_stats->SetX2NDC(0.9);
    fake_rate_project_y_stats->SetY2NDC(0.99);
  }
  pad4->SetLeftMargin(0.15);
  pad4->SetRightMargin(0.1);
  pad4->SetBottomMargin(0.17);
  pad4->SetTopMargin(0.1);
  pad4->Update();

  TString picture_name = "fake_rate_";
  picture_name.Append(numerator);
  picture_name.Append("_");
  picture_name.Append(sample);
  picture_name.Append(".gif");
  picture_name.ReplaceAll("num_","");

  canvas.Print(picture_name);

}

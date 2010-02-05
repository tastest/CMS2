//void overlay4(TH1F* h1, TH1F* h2, TH1F* h3, TH1F* h4, char* title){
void overlay2(TH1F* h1, TH1F* h2, const char* plotname){

  gStyle->SetOptStat("nemroui");
  //	string outfile = title;
  //	outfile.replace(3, 1, "_");
  //	outfile.replace( outfile.find(")"), 2, "_");
  //	outfile.replace( outfile.find("-") - 1, 3, "_");
  //	outfile.replace( outfile.find("Selection") - 1, outfile.length(), ".jpg");
  //	string outfile_log = outfile;
  //	outfile_log.replace( outfile.find(".jpg"), outfile.length(), "_log.jpg");
  //
  // get maximum range
  double ymax = h1->GetMaximum();
  if(h2->GetMaximum() > ymax){
    ymax = h2->GetMaximum();
  }
  //	if(h3->GetMaximum() > ymax){
  //		ymax = h3->GetMaximum();
  //	}
  //	if(h4->GetMaximum() > ymax){
  //		ymax = h4->GetMaximum();
  //	}
  //
  // 	// draw histograms
  gStyle->SetTitleW( .57 );
  TString buffer = "c1_";
  buffer += plotname;
  TCanvas *c1 = new TCanvas(buffer,buffer,1280,960);
  //	c1->SetTitle(title);
  //	h1->SetTitle(title);

  h1->SetMarkerStyle(3);
  h2->SetMarkerStyle(25);


  h1->Draw();
  h2->Draw("sames");
  //  h3->Draw("sames");
  //  h4->Draw("sames");
  //
  //	// format style
  h1->SetMaximum(ymax*1.1);
  //  h1->SetXTitle("MT2");
  //  h2->SetXTitle("MT2");
  //  h3->SetXTitle("MT2");
  //  h4->SetXTitle("MT2");
  gPad->Update();
  //h1->SetLineColor(kBlack);
  //h2->SetLineColor(kRed);
  //	h3->SetLineColor(kRed);
  //	h4->SetLineColor(kGreen);
  TPaveStats* st1 = (TPaveStats*)h1->FindObject("stats");
  TPaveStats* st2 = (TPaveStats*)h2->FindObject("stats");
  //	TPaveStats* st3 = (TPaveStats*)h3->FindObject("stats");
  //	TPaveStats* st4 = (TPaveStats*)h4->FindObject("stats");
  st1->SetTextColor(kBlack);
  st2->SetTextColor(kRed);
  //	st3->SetTextColor(kRed);
  //	st4->SetTextColor(kGreen);
  st2->SetX1NDC( 2*st1->GetX1NDC() - st1->GetX2NDC() );
  st2->SetX2NDC( st1->GetX1NDC() );
  //	st3->SetY1NDC( 2*st1->GetY1NDC() - st1->GetY2NDC() );
  //	st3->SetY2NDC( st1->GetY1NDC() );
  //	st4->SetX1NDC( 2*st1->GetX1NDC() - st1->GetX2NDC() );
  //	st4->SetX2NDC( st1->GetX1NDC() );
  //	st4->SetY1NDC( 2*st1->GetY1NDC() - st1->GetY2NDC() );
  //	st4->SetY2NDC( st1->GetY1NDC() );
  st1->Draw();
  st2->Draw();
  //	st3->Draw();
  //	st4->Draw();
  //	h1->SetName("Modified All Jets no #eta cut");
  //	h2->SetName("Leading 2 Jets");
  //	h3->SetName("All Jets");
  //	h4->SetName("All Jets no #eta cut");
  gPad->Update();
  //	c1->SaveAs( Form("../plots/jpg/%s",outfile.c_str()) );
  //	c1->SetLogy();
  //	c1->SaveAs( Form("../plots/jpg/%s",outfile_log.c_str()) );
  c1->SaveAs( Form("%s",plotname) );
  return;
}

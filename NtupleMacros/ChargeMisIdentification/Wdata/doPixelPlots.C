{

  TCanvas *canvas = new TCanvas;
  canvas->Divide(3,2);

  TPad *pad1 = (TPad*)canvas->cd(1);
  pad1->SetLeftMargin(0.16);
  pad1->SetBottomMargin(0.16);

  els_chargeFirstPixelHit_corCharge->Draw("HIST");

  TPad *pad4 = (TPad*)canvas->cd(4);
  pad4->SetLeftMargin(0.16);
  pad4->SetBottomMargin(0.16);

  els_chargeFirstPixelHit_incorCharge->Draw("HIST");

  TPad *pad2 = (TPad*)canvas->cd(2);
  pad2->SetLeftMargin(0.16);
  pad2->SetBottomMargin(0.16);

  els_chargeFirstPixelHit_corCharge_barrel->Draw("HIST");

  TPad *pad5 = (TPad*)canvas->cd(5);
  pad5->SetLeftMargin(0.16);
  pad5->SetBottomMargin(0.16);

  els_chargeFirstPixelHit_incorCharge_barrel->Draw("HIST");

  TPad *pad3 = (TPad*)canvas->cd(3);
  pad3->SetLeftMargin(0.16);
  pad3->SetBottomMargin(0.16);

  els_chargeFirstPixelHit_corCharge_forward->Draw("HIST");

  TPad *pad6 = (TPad*)canvas->cd(6);
  pad6->SetLeftMargin(0.16);
  pad6->SetBottomMargin(0.16);

  els_chargeFirstPixelHit_incorCharge_forward->Draw("HIST");

}

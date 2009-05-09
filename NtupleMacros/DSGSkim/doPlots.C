{

  TCanvas *canvas = new TCanvas;

  TH1F *els_pt = (TH1F*)gFile->Get("els_pt");
  TH1F *els_eta = (TH1F*)gFile->Get("els_eta");
  TH1F *els_pt_corCharge = (TH1F*)gFile->Get("els_pt_corCharge");
  TH1F *els_eta_corCharge = (TH1F*)gFile->Get("els_eta_corCharge");
  TH1F *els_pt_incorCharge = (TH1F*)gFile->Get("els_pt_incorCharge");
  TH1F *els_eta_incorCharge = (TH1F*)gFile->Get("els_eta_incorCharge");

  canvas->Clear();
  els_pt->Draw();
  canvas->Print("els_pt.png");

  canvas->Clear();
  els_eta->Draw();
  canvas->Print("els_eta.png");

  canvas->Clear();
  els_pt_incorCharge->Divide(els_pt);
  els_pt_incorCharge->Draw();
  canvas->Print("els_pt_wrongChargeFraction.png");

  canvas->Clear();
  els_eta_incorCharge->Divide(els_eta);
  els_eta_incorCharge->Draw();
  canvas->Print("els_eta_wrongChargeFraction.png");


}

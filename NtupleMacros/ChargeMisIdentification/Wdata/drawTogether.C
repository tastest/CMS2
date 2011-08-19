void drawTogether()
{
//charflipRej(){
gStyle->SetOptStat(0);
TFile* file_336p4 = TFile::Open("efficiency_3x.root");

TFile* file_361p3 = TFile::Open("efficiency_126_Jura_nohitpattern.root");
TCanvas *canvas = new TCanvas;
TPaveText *box = 0;

canvas->Clear();
canvas->SetLeftMargin(0.155);
canvas->SetBottomMargin(0.155);
TH1* reco_eff_pt_336p4  = file_336p4->Get("RecoEff_pt"); 
TH1* reco_eff_pt_361p3  = file_361p3->Get("RecoEff_pt"); 
reco_eff_pt_361p3->SetLineColor(kRed);
reco_eff_pt_361p3->SetMarkerColor(kRed);
reco_eff_pt_336p4->Draw();
reco_eff_pt_361p3->Draw("same");
canvas->Print("RecoEfficiencyPt_336p4_361p3.png");

box = new TPaveText(40.,0.4,80.,0.5);
box->SetBorderSize(1);
box->SetFillColor(0);
box->InsertText(Form("Red: 361p3; Black: 336p4"));
box->Draw();
canvas->Print("RecoEfficiencyPt_336p4_361p3.png");

canvas->Clear();
canvas->SetLeftMargin(0.155);
canvas->SetBottomMargin(0.155);
TH1* reco_eff_eta_336p4  = file_336p4->Get("RecoEff_eta"); 
TH1* reco_eff_eta_361p3  = file_361p3->Get("RecoEff_eta"); 
reco_eff_eta_361p3->SetLineColor(kRed);
reco_eff_eta_361p3->SetMarkerColor(kRed);
reco_eff_eta_336p4->Draw();
reco_eff_eta_361p3->Draw("same");
canvas->Print("RecoEfficiencyEta_336p4_361p3.png");

box = new TPaveText(40.,0.4,80.,0.5);
box->SetBorderSize(1);
box->SetFillColor(0);
box->InsertText(Form("Red: 361p3; Black: 336p4"));
box->Draw();
canvas->Print("RecoEfficiencyEta_336p4_361p3.png");

canvas->Clear();
canvas->SetLeftMargin(0.155);
canvas->SetBottomMargin(0.155);
TH1* ChargeMisIdRate_eta_336p4  = file_336p4->Get("ChargeMisIdRate_eta"); 
TH1* ChargeMisIdRate_eta_361p3  = file_361p3->Get("ChargeMisIdRate_eta"); 
ChargeMisIdRate_eta_361p3->SetLineColor(kRed);
ChargeMisIdRate_eta_361p3->SetMarkerColor(kRed);
ChargeMisIdRate_eta_336p4->Draw();
ChargeMisIdRate_eta_361p3->Draw("same");
canvas->Print("ChargeMisIdRateEta_336p4_361p3.png");

box = new TPaveText(40.,0.4,80.,0.5);
box->SetBorderSize(1);
box->SetFillColor(0);
box->InsertText(Form("Red: 361p3; Black: 336p4"));
box->Draw();
canvas->Print("ChargeMisIdRateEta_336p4_361p3.png");

canvas->Clear();
canvas->SetLeftMargin(0.155);                     
canvas->SetBottomMargin(0.155);
TH1* ChargeMisIdRate_pt_336p4  = (TH1*)file_336p4->Get("ChargeMisIDRate_pt"); 
TH1* ChargeMisIdRate_pt_361p3  = (TH1*)file_361p3->Get("ChargeMisIDRate_pt"); 
ChargeMisIdRate_pt_361p3->SetLineColor(kRed);
ChargeMisIdRate_pt_361p3->SetMarkerColor(kRed);
ChargeMisIdRate_pt_336p4->Draw();
ChargeMisIdRate_pt_361p3->Draw("same");
canvas->Print("ChargeMisIdRatePt_336p4_361p3.png");

box = new TPaveText(40.,0.4,80.,0.5);
box->SetBorderSize(1);
box->SetFillColor(0);
box->InsertText(Form("Red: 361p3; Black: 336p4"));
box->Draw();
canvas->Print("ChargeMisIdRatePt_336p4_361p3.png");

}

/*
after_charflipRej(){
gStyle->SetOptStat(0);
TFile* file_336p4 = TFile::Open("/home/users/yanjuntu/CMSSW_3_3_6/src/ChargeFlip_336p4/CMS2/NtupleMacros/ChargeMisIdentification/efficiency_336p4_after_charflipRej.root");

TFile* file_361p3 = TFile::Open("/home/users/yanjuntu/CMSSW_3_3_6/src/ChargeFlip/CMS2/NtupleMacros/ChargeFlip/Ana_looper/efficiency_361p3_after_charflipRej_can02.root");
TCanvas *canvas = new TCanvas;
TPaveText *box = 0;

canvas->Clear();
canvas->SetLeftMargin(0.155);
canvas->SetBottomMargin(0.155);
TH1* reco_eff_pt_336p4  = file_336p4->Get("RecoEff_pt"); 
TH1* reco_eff_pt_361p3  = file_361p3->Get("RecoEff_pt"); 
reco_eff_pt_361p3->SetLineColor(kRed);
reco_eff_pt_361p3->SetMarkerColor(kRed);
reco_eff_pt_336p4->Draw();
reco_eff_pt_361p3->Draw("same");
canvas->Print("RecoEfficiencyPt_336p4_361p3_after.png");

box = new TPaveText(40.,0.4,80.,0.5);
box->SetBorderSize(1);
box->SetFillColor(0);
box->InsertText(Form("Red: 361p3; Black: 336p4"));
box->Draw();
canvas->Print("RecoEfficiencyPt_336p4_361p3_after.png");

canvas->Clear();
canvas->SetLeftMargin(0.155);
canvas->SetBottomMargin(0.155);
TH1* reco_eff_eta_336p4  = file_336p4->Get("RecoEff_eta"); 
TH1* reco_eff_eta_361p3  = file_361p3->Get("RecoEff_eta"); 
reco_eff_eta_361p3->SetLineColor(kRed);
reco_eff_eta_361p3->SetMarkerColor(kRed);
reco_eff_eta_336p4->Draw();
reco_eff_eta_361p3->Draw("same");
canvas->Print("RecoEfficiencyEta_336p4_361p3_after.png");

box = new TPaveText(40.,0.4,80.,0.5);
box->SetBorderSize(1);
box->SetFillColor(0);
box->InsertText(Form("Red: 361p3; Black: 336p4"));
box->Draw();
canvas->Print("RecoEfficiencyEta_336p4_361p3_after.png");

canvas->Clear();
canvas->SetLeftMargin(0.155);
canvas->SetBottomMargin(0.155);
TH1* ChargeMisIdRate_eta_336p4  = file_336p4->Get("ChargeMisIdRate_eta"); 
TH1* ChargeMisIdRate_eta_361p3  = file_361p3->Get("ChargeMisIdRate_eta"); 
ChargeMisIdRate_eta_361p3->SetLineColor(kRed);
ChargeMisIdRate_eta_361p3->SetMarkerColor(kRed);
ChargeMisIdRate_eta_336p4->Draw();
ChargeMisIdRate_eta_361p3->Draw("same");
canvas->Print("ChargeMisIdRateEta_336p4_361p3_after.png");

box = new TPaveText(40.,0.4,80.,0.5);
box->SetBorderSize(1);
box->SetFillColor(0);
box->InsertText(Form("Red: 361p3; Black: 336p4"));
box->Draw();
canvas->Print("ChargeMisIdRateEta_336p4_361p3_after.png");


canvas->Clear();
canvas->SetLeftMargin(0.155);                     
canvas->SetBottomMargin(0.155);
TH1* ChargeMisIdRate_pt_336p4  = (TH1*)file_336p4->Get("ChargeMisIDRate_pt"); 
TH1* ChargeMisIdRate_pt_361p3  = (TH1*)file_361p3->Get("ChargeMisIDRate_pt"); 
ChargeMisIdRate_pt_361p3->SetLineColor(kRed);
ChargeMisIdRate_pt_361p3->SetMarkerColor(kRed);
ChargeMisIdRate_pt_336p4->Draw();
ChargeMisIdRate_pt_361p3->Draw("same");
canvas->Print("ChargeMisIdRatePt_336p4_361p3_after.png");

box = new TPaveText(40.,0.4,80.,0.5);
box->SetBorderSize(1);
box->SetFillColor(0);
box->InsertText(Form("Red: 361p3; Black: 336p4"));
box->Draw();
canvas->Print("ChargeMisIdRatePt_336p4_361p3_after.png");

}
*/

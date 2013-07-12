/* 
   Code to do the data driven DY normalization for you.
   Takes in 4 mandatory variables:
   
   dataFName = path with histograms (NJet mainly)
   dataNoMETFName = path with histograms with NO MET cut. This is crucial
   mcFName = path of MC file that has the same cuts applied as dataFName
   scaleMC = lumi to scale MC to
   
   verbose = provides dump of variables used for DY estimate such as
   R_out/in, k_ee, k_mm, etc
*/


#include "TH1F.h"
#include "TString.h"
#include "histtools.C"
#include <iostream>
#include "TFile.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "CommonFunctions.C"
using namespace std;

TH1F* DYee_hnJetoutZwindow_ee;
TH1F* DYee_hnJetinZwindow_ee;
TH1F* data_hnJetinZwindow_ee;
TH1F* data_hnJetoutZwindow_ee;

TH1F* DYmm_hnJetoutZwindow_mm;
TH1F* DYmm_hnJetinZwindow_mm;
TH1F* data_hnJetinZwindow_mm;
TH1F* data_hnJetoutZwindow_mm;

TH1F* data_hnJetinZwindow_em;

TH1F* data_hnJetinZwindow_ee_noMET;
TH1F *data_hnJetinZwindow_mm_noMET;

void getEstimates(int lowJetBin, int highJetBin, 
		  float& n_ee, float& n_mm, 
		  float& n_eeE, float& n_mmE,
		  bool verbose, float scaleMC,
		  char* formatS = "%6.3f") {

  
  //get the loose MET cut hists    
  float MCo_ee	= GetEntries(DYee_hnJetoutZwindow_ee, lowJetBin, highJetBin);
  float MCi_ee	= GetEntries(DYee_hnJetinZwindow_ee, lowJetBin, highJetBin);
  float MCo_eeE	= GetTotalError(DYee_hnJetoutZwindow_ee, lowJetBin, highJetBin);
  float MCi_eeE	= GetTotalError(DYee_hnJetinZwindow_ee, lowJetBin, highJetBin);
  float Do_ee	= GetEntries(data_hnJetoutZwindow_ee, lowJetBin, highJetBin);
  float Di_ee	= GetEntries(data_hnJetinZwindow_ee, lowJetBin, highJetBin);
  float Di_eeE	= sqrt(Di_ee);
  float R_ee	= MCo_ee/MCi_ee;
  float R_eeE	= sqrt(pow(MCo_eeE/MCi_ee,2) + pow(MCo_ee*MCi_eeE/MCi_ee/MCi_ee, 2)/4.);
    

  float MCo_mm	= GetEntries(DYmm_hnJetoutZwindow_mm, lowJetBin, highJetBin);
  float MCi_mm	= GetEntries(DYmm_hnJetinZwindow_mm, lowJetBin, highJetBin);
  float MCo_mmE	= GetTotalError(DYmm_hnJetoutZwindow_mm, lowJetBin, highJetBin);
  float MCi_mmE	= GetTotalError(DYmm_hnJetoutZwindow_mm, lowJetBin, highJetBin);
  float Do_mm	= GetEntries(data_hnJetoutZwindow_mm, lowJetBin, highJetBin);
  float Di_mm	= GetEntries(data_hnJetinZwindow_mm, lowJetBin, highJetBin);
  float Di_mmE	= sqrt(Di_mm);
  float R_mm	= MCo_mm/MCi_mm;
  float R_mmE	= sqrt(pow(MCo_mmE/MCi_mm,2) + pow(MCo_mm*MCi_mmE/MCi_mm/MCi_mm, 2)/4.);


  //Nentries, emu
  float Di_em		= GetEntries(data_hnJetinZwindow_em, lowJetBin, highJetBin);
  float Di_emE	= sqrt(Di_em);
    
  //Nentries, no MET cut
  float Di_ee_noMET	= GetEntries(data_hnJetinZwindow_ee_noMET, lowJetBin, highJetBin);
  float Di_eeE_noMET	= sqrt(Di_ee_noMET);
  float Di_mm_noMET	= GetEntries(data_hnJetinZwindow_mm_noMET, lowJetBin, highJetBin);
  float Di_mmE_noMET	= sqrt(Di_mm_noMET);
  
  float k_ee		= sqrt(Di_ee_noMET/Di_mm_noMET);
  float k_eeE		= sqrt(pow(Di_eeE_noMET/Di_mm_noMET,2) + pow(Di_mmE_noMET*Di_ee_noMET/Di_mm_noMET/Di_mm_noMET,2)/4.);
  
  float k_mm		= 1/k_ee;
  float k_mmE		= sqrt(pow(Di_mmE_noMET/Di_ee_noMET,2) + pow(Di_eeE_noMET*Di_mm_noMET/Di_ee_noMET/Di_ee_noMET,2)/4.);
  

  float pred_ee	= R_ee*(Di_ee - 0.5*Di_em*k_ee);
  float pred_eeE	= sqrt(pow((Di_ee-0.5*Di_em*k_ee)*R_eeE,2) + pow(R_ee*Di_eeE,2) + pow(0.5*R_ee*k_ee*Di_emE,2) + pow(0.5*R_ee*Di_em*k_eeE,2));

  float pred_mm	= R_mm*(Di_mm - 0.5*Di_em*k_mm);
  float pred_mmE	= sqrt(pow((Di_mm-0.5*Di_em*k_mm)*R_mmE,2) + pow(R_mm*Di_mmE,2) + pow(0.5*R_mm*k_mm*Di_emE,2) + pow(0.5*R_mm*Di_em*k_mmE,2));


  if(verbose) {
    cout << "Di_ee: " << Di_ee << "+/-" << Di_eeE << endl;
    cout << "R_ee: " << R_ee << "+/-" << R_eeE << endl;
    cout << "k_ee: " << k_ee << "+/-" << k_eeE << endl;
    cout << "Di_em: " << Di_em << "+/-" << Di_emE << endl;

    /*
      cout << "Prediction from MC in ee channel: " << scaleMC*GetEntries(DYee_hnJetoutZwindow_ee, lowJetBin, highJetBin) 
      << " \\pm " << scaleMC*GetTotalError(DYee_hnJetoutZwindow_ee, lowJetBin, highJetBin) << endl;
      if(MCi_ee > 0)
      cout << "Number predicted in Data in ee: " << pred_ee << " \\pm " << pred_ee << endl;
      else 
      cout << "Number predicted in Data in ee: " << 0 << endl;
      cout << "Number of events in Data outside of the Z region in the ee channel: " << Do_ee << endl;
      cout << endl;
    */

    cout << "Di_mm: " << Di_mm << "+/-" << Di_mmE << endl;
    cout << "R_mm: " << R_mm << "+/-" << R_mmE << endl;
    cout << "k_mm: " << k_mm << "+/-" << k_mmE << endl;
    cout << "Di_em: " << Di_em << "+/-" << Di_emE << endl;


    cout << endl << endl;
    cout << "\\begin{tabular}{l |  c  c  c  c}" << endl; 
    cout << "\\hline" << endl;
    cout << "Simple & ee & $\\mu\\mu$ & e$\\mu$ & all \\\\ \\hline" << endl;
    cout << "MC Prediction & " << formatFloat(scaleMC*GetEntries(DYee_hnJetoutZwindow_ee, lowJetBin, highJetBin), formatS)
	 << " $\\pm$ " << formatFloat(scaleMC*GetTotalError(DYee_hnJetoutZwindow_ee, lowJetBin, highJetBin), formatS) << " & ";
    cout << formatFloat(scaleMC*GetEntries(DYmm_hnJetoutZwindow_mm, lowJetBin, highJetBin),formatS) 
	 << " $\\pm$ " << formatFloat(scaleMC*GetTotalError(DYmm_hnJetoutZwindow_mm, lowJetBin, highJetBin),formatS) << " & - & ";
    float blah = scaleMC*(GetEntries(DYee_hnJetoutZwindow_ee, lowJetBin,highJetBin) + 
			  GetEntries(DYmm_hnJetoutZwindow_mm, lowJetBin, highJetBin));
    float blahE = scaleMC*(sqrt(pow(GetTotalError(DYee_hnJetoutZwindow_ee, lowJetBin,highJetBin), 2) + 
				pow(GetTotalError(DYmm_hnJetoutZwindow_mm, lowJetBin, highJetBin), 2)));
    cout << formatFloat(blah,formatS) << "$\\pm$ " << formatFloat(blahE,formatS) << " \\\\ " << endl;
    if(MCi_ee > 0)
      cout << "Data Driven DY Estimate & " << formatFloat(pred_ee,formatS) << " $\\pm$ " << formatFloat(pred_eeE,formatS) << " & ";
    else 
      cout << "Data Driven DY Estimate & " << 0 << " $\\pm$ " << 0 << " & ";
    if(MCi_mm > 0)
      cout << formatFloat(pred_mm,formatS) << " $\\pm$ " << formatFloat(pred_mmE, formatS);
    else 
      cout << 0 << " $\\pm$ " << 0;//<< "\\\\ \\hline" << endl;
    cout << " & - & " << formatFloat(pred_ee + pred_mm, formatS) << " $\\pm$ " << formatFloat(sqrt(pred_eeE*pred_eeE + pred_mmE*pred_mmE),formatS) << "\\\\ \\hline" << endl;
    cout << "Actual Yield in Data & " << Do_ee << " & " << Do_mm << " & - & " << Do_ee + Do_mm << "\\\\ \\hline" << endl;
    cout << "\\end{tabular}" << endl;


    /*
      cout << "Prediction from MC in mm channel: " << scaleMC*GetEntries(DYmm_hnJetoutZwindow_mm, lowJetBin, highJetBin) 
      << " \\pm " << scaleMC*GetTotalError(DYmm_hnJetoutZwindow_mm, lowJetBin, highJetBin) << endl;  
      if(MCi_mm > 0)
      cout << "Number predicted in Data in mm: " << pred_mm << " \\pm " << pred_mmE << endl;
      else 
      cout << "Number predicted in Data in mm: " << 0 << endl;
    */

  }
  n_ee = pred_ee;
  n_eeE = pred_eeE;
  n_mm = pred_mm;
  n_mmE = pred_mmE;
    
}





void DYest(TString dataFName, TString dataNoMETFName, TString mcFName, float scaleMC, bool verbose = true){

  TFile *f_data		= TFile::Open(dataFName.Data());
  TFile *f_mc		= TFile::Open(mcFName.Data());
  TFile *f_dataNoMET	= TFile::Open(dataNoMETFName.Data());
  
  
  DYee_hnJetoutZwindow_ee	= (TH1F*)f_mc->Get("DYee_hnJetoutZwindow_ee");
  DYee_hnJetinZwindow_ee	= (TH1F*)f_mc->Get("DYee_hnJetinZwindow_ee");
  data_hnJetinZwindow_ee	= (TH1F*)f_data->Get("data_hnJetinZwindow_ee");
  data_hnJetoutZwindow_ee	= (TH1F*)f_data->Get("data_hnJetoutZwindow_ee");
  
  DYmm_hnJetoutZwindow_mm	= (TH1F*)f_mc->Get("DYmm_hnJetoutZwindow_mm");
  DYmm_hnJetinZwindow_mm	= (TH1F*)f_mc->Get("DYmm_hnJetinZwindow_mm");
  data_hnJetinZwindow_mm	= (TH1F*)f_data->Get("data_hnJetinZwindow_mm");
  data_hnJetoutZwindow_mm	= (TH1F*)f_data->Get("data_hnJetoutZwindow_mm");


  //get the emu data hists
  data_hnJetinZwindow_em	= (TH1F*)f_data->Get("data_hnJetinZwindow_em");

  //get the histos from the loose cuts
  data_hnJetinZwindow_ee_noMET	= (TH1F*)f_dataNoMET->Get("data_hnJetinZwindow_ee");
  data_hnJetinZwindow_mm_noMET	= (TH1F*)f_dataNoMET->Get("data_hnJetinZwindow_mm");
  
  //float temp_ee = 0;
  //float temp_mm = 0;

  float n_ee = 0;
  float n_eeE =0; 
  float n_mm = 0;
  float n_mmE =0;


  const char *jetbins[5] = {"0", "1", "2", "3", "#geq 4"};
  TH1F *DYEst_hnJet_ee = new TH1F("DYEst_hnJet_ee", "Data Driven DY->ee estimate", 
				  DYee_hnJetinZwindow_ee->GetNbinsX(), 
				  DYee_hnJetinZwindow_ee->GetXaxis()->GetXmin(),
				  DYee_hnJetinZwindow_ee->GetXaxis()->GetXmax());
  TH1F *DYEst_hnJet_mm = new TH1F("DYEst_hnJet_mm", "Data Driven DY->mm estimate", 
				  DYee_hnJetinZwindow_ee->GetNbinsX(), 
				  DYee_hnJetinZwindow_ee->GetXaxis()->GetXmin(),
				  DYee_hnJetinZwindow_ee->GetXaxis()->GetXmax());
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  DYEst_hnJet_ee->Sumw2();
  DYEst_hnJet_mm->Sumw2();
  DYEst_hnJet_ee->SetDirectory(rootdir);
  DYEst_hnJet_mm->SetDirectory(rootdir);


  for(int k = 0; k<5; k++) {
    DYEst_hnJet_ee->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
    DYEst_hnJet_ee->GetXaxis()->SetLabelSize(0.07);

    DYEst_hnJet_mm->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
    DYEst_hnJet_mm->GetXaxis()->SetLabelSize(0.07);
  }
  
  cout << "\\documentclass{article}" << endl;
  cout << "\\usepackage{times}" << endl;
  cout << "\\begin{document}" << endl;


  if(verbose) {
    cout << endl << "//-----------------------------------------------------------------------//" << endl << endl;
    cout << "Doing the 0 jet Bin:" << endl;
  }
  getEstimates(1, 1, n_ee, n_mm, n_eeE, n_mmE, verbose, scaleMC);
  DYEst_hnJet_ee->SetBinContent(1, n_ee);
  DYEst_hnJet_mm->SetBinContent(1, n_mm);
  DYEst_hnJet_ee->SetBinError(1, n_eeE);
  DYEst_hnJet_mm->SetBinError(1, n_mmE);
  //temp_ee = temp_ee + n_ee;
  //temp_mm = temp_mm + n_mm;
  if(verbose) {
    cout << endl << "//-----------------------------------------------------------------------//" << endl << endl;
    cout << "Doing the 1 jet Bin: " << endl;
  }
  getEstimates(2, 2, n_ee, n_mm, n_eeE, n_mmE, verbose, scaleMC);
  DYEst_hnJet_ee->SetBinContent(2, n_ee);
  DYEst_hnJet_mm->SetBinContent(2, n_mm);
  DYEst_hnJet_ee->SetBinError(2, n_eeE);
  DYEst_hnJet_mm->SetBinError(2, n_mmE);
  //temp_ee = temp_ee + n_ee;
  //temp_mm = temp_mm + n_mm;
  if(verbose) {
    cout << endl << "//-----------------------------------------------------------------------//" << endl << endl;
    cout << "Doing the >= 2 jet Bin: " << endl;
  }
  getEstimates(3, 5, n_ee, n_mm, n_eeE, n_mmE, verbose, scaleMC);

  /*
  if(verbose) {
    cout << endl << "//-----------------------------------------------------------------------//" << endl << endl;
    cout << "Doing the 2 jet Bin: " << endl;
  }
  */
  getEstimates(3, 3, n_ee, n_mm, n_eeE, n_mmE, false, scaleMC);

  DYEst_hnJet_ee->SetBinContent(3, n_ee);
  DYEst_hnJet_mm->SetBinContent(3, n_mm);
  DYEst_hnJet_ee->SetBinError(3, n_eeE);
  DYEst_hnJet_mm->SetBinError(3, n_mmE);
  //temp_ee = temp_ee + n_ee;
  //temp_mm = temp_mm + n_mm;

  /*
  if(verbose) {
    cout << endl << "//-----------------------------------------------------------------------//" << endl << endl;
    cout << "Doing the 3 jet Bin: " << endl;
  }
  */
  getEstimates(4, 4, n_ee, n_mm, n_eeE, n_mmE, false, scaleMC);

  DYEst_hnJet_ee->SetBinContent(4, n_ee);
  DYEst_hnJet_mm->SetBinContent(4, n_mm);
  DYEst_hnJet_ee->SetBinError(4, n_eeE);
  DYEst_hnJet_mm->SetBinError(4, n_mmE);
  //temp_ee = temp_ee + n_ee;
  //temp_mm = temp_mm + n_mm;

  /*
  if(verbose) {
    cout << endl << "//-----------------------------------------------------------------------//" << endl << endl;
    cout << "Doing the >= 4 jet Bin: " << endl;
  }
  */
  getEstimates(5, 5, n_ee, n_mm, n_eeE, n_mmE, false, scaleMC);

  DYEst_hnJet_ee->SetBinContent(5, n_ee);
  DYEst_hnJet_mm->SetBinContent(5, n_mm);
  DYEst_hnJet_ee->SetBinError(5, n_eeE);
  DYEst_hnJet_mm->SetBinError(5, n_mmE);
  //temp_ee = temp_ee + n_ee;
  //temp_mm = temp_mm + n_mm;
  if(verbose) {
    cout << endl << "//-----------------------------------------------------------------------//" << endl << endl;
    //cout << "total ee: " << temp_ee << endl;
    //cout << "total mm: " << temp_mm << endl;
    cout << "Inclusive DY estimate: " << endl;
  }
  getEstimates(1, 5, n_ee, n_mm, n_eeE, n_mmE, verbose, scaleMC);

  cout << "\\end{document}" << endl;
    
  //hist::deleteHistos("DYee_*");

  f_data->Close();
  f_mc->Close();
  f_dataNoMET->Close();
  
}


  

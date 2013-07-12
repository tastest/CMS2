#include <iostream>
#include <utility>
#include "TFile.h"
#include "TH1F.h"
#include "TSystem.h"

#include "CommonFunctions.C"

//dataFName - the data file with the set of cuts that you want the FRs to be estimated for
//dataFName is NOT the Fake Rate file!!!
//The QCD and WJets FR files are assumed using the dataFName
//useZmassRegionForSpillage - Do you want to estimate the spillage using the Zmass region? 
//With low stat, this only makes sense if you have not applied the MET or 2 jets requirements yet. 


vector<TH1F*> EstimateSignalSpillage(TString dataFName, bool combineHigherJetMultBins, bool useZmassRegionForSpillage, bool verbose = true, const char* formatS = "%6.4f") { 
	
  if(combineHigherJetMultBins) {
    cout << "//---------------------------------------------------------------------------//" << endl << endl;
    cout << "You are choosing to get results after combining higher jet multiplicity bins. " << endl;
    cout << "The Spillage Fraction histograms will not make sense, so don't use them!!!!" << endl;
    cout << "//---------------------------------------------------------------------------//" << endl << endl;
  }

  if(dataFName.Contains("FRhist")) {
    cout << "The first argument should not be the one with the Fake Rate yields, but the full yields after the cuts you are investigating" << endl;
    cout << "The code will figure out the FR files from dataFName. " << endl;
    cout << "Exiting" << endl;
  }
  
  bool haveAppliedZvetoCut = false;
  if(dataFName.Contains("vetoZmass") )
    haveAppliedZvetoCut = true;


  //get the FR files
  TString wJetsFile = dataFName;
  TString qcdFile = dataFName;
  wJetsFile.ReplaceAll("hist", "FRhist");
  wJetsFile.ReplaceAll(".root", "_estimateWJets.root");
  wJetsFile.ReplaceAll("applylepIDCuts_applylepIsoCuts", "applyFOv2Cuts");
  
  
  qcdFile.ReplaceAll("hist", "FRhist");
  qcdFile.ReplaceAll(".root", "_estimateQCD.root");
  qcdFile.ReplaceAll("applylepIDCuts_applylepIsoCuts", "applyFOv2Cuts");

  //return convention is bizarre
  if(gSystem->AccessPathName(wJetsFile.Data()) == 1) {
    cout << "Cannot find WJets FR file: " << wJetsFile << endl;
    cout << "Exiting" << endl;
    gSystem->Exit(1);
  }
  if(gSystem->AccessPathName(qcdFile.Data()) == 1) {
    cout << "Cannot find QCD FR file: " << qcdFile << endl;
    cout << "Exiting" << endl;
    gSystem->Exit(1);
  }

  vector<TH1F*> v_SR;
  v_SR.reserve(3);

  TFile *f_Data    = TFile::Open(dataFName.Data(), "READ");
  TFile *f_FRwJets = TFile::Open(wJetsFile.Data(), "READ");
  TFile *f_FRqcd = TFile::Open(qcdFile.Data(), "READ");
  
  TH1F *h_data_ee;
  TH1F *h_data_mm;
  TH1F *h_data_em = (TH1F*)f_Data->Get("data_hnJet_em");

  TH1F *h_wJets_FRee;
  TH1F *h_wJets_FRmm;
  TH1F *h_wJets_FRem = (TH1F*)f_FRwJets->Get("data_hnJet_em");  
  
  TH1F *h_qcd_FRee;
  TH1F *h_qcd_FRmm;
  TH1F *h_qcd_FRem = (TH1F*)f_FRqcd->Get("data_hnJet_em");  
  
  if(useZmassRegionForSpillage) {
    h_data_ee = (TH1F*)f_Data->Get("data_hnJetinZwindow_ee");
    h_data_mm = (TH1F*)f_Data->Get("data_hnJetinZwindow_mm");
    h_wJets_FRee = (TH1F*)f_FRwJets->Get("data_hnJetinZwindow_ee");
    h_wJets_FRmm = (TH1F*)f_FRwJets->Get("data_hnJetinZwindow_mm");
    h_qcd_FRee = (TH1F*)f_FRqcd->Get("data_hnJetinZwindow_ee");
    h_qcd_FRmm = (TH1F*)f_FRqcd->Get("data_hnJetinZwindow_mm");
  } else {
    h_data_ee = (TH1F*)f_Data->Get("data_hnJet_ee");
    h_data_mm = (TH1F*)f_Data->Get("data_hnJet_mm");
    h_wJets_FRee = (TH1F*)f_FRwJets->Get("data_hnJet_ee");
    h_wJets_FRmm = (TH1F*)f_FRwJets->Get("data_hnJet_mm");
    h_qcd_FRee = (TH1F*)f_FRqcd->Get("data_hnJet_ee");
    h_qcd_FRmm = (TH1F*)f_FRqcd->Get("data_hnJet_mm");
  }
  


  TH1F *h_FR_ee = (TH1F*)h_wJets_FRee->Clone();
  TH1F *h_FR_mm = (TH1F*)h_wJets_FRmm->Clone();
  TH1F *h_FR_em = (TH1F*)h_wJets_FRem->Clone();

 
  if(h_FR_ee->GetSumw2()->GetSize() == 0)
    h_FR_ee->Sumw2();
  if(h_FR_mm->GetSumw2()->GetSize() == 0)
    h_FR_mm->Sumw2();
  if(h_FR_em->GetSumw2()->GetSize() == 0)
    h_FR_em->Sumw2();

  
  TH1F *h_SR_ee = (TH1F*)h_FR_ee->Clone();
  TH1F *h_SR_mm = (TH1F*)h_FR_mm->Clone();
  TH1F *h_SR_em = (TH1F*)h_FR_em->Clone();

  h_SR_ee->SetName("h_SR_ee");
  h_SR_mm->SetName("h_SR_mm");
  h_SR_em->SetName("h_SR_em");

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  h_SR_ee->SetDirectory(rootdir);
  h_SR_mm->SetDirectory(rootdir);
  h_SR_em->SetDirectory(rootdir);

  string colSep = " & ";
  string pm = " $\\pm$ ";
  string endL = " \\\\";

  
  double nWJets_ee_total = 0;
  double nWJets_eeE_total = 0;
  double nQCD_ee_total = 0;
  double nQCD_eeE_total = 0;

  double nWJets_mm_total = 0;
  double nWJets_mmE_total = 0;
  double nQCD_mm_total = 0;
  double nQCD_mmE_total = 0;

  double nWJets_em_total = 0;
  double nWJets_emE_total = 0;
  double nQCD_em_total = 0;
  double nQCD_emE_total = 0;


  
  for(int ibinx = 1; ibinx < 6; ibinx++) {

    int lowBin = ibinx;
    int highBin = ibinx;
    
    if(combineHigherJetMultBins) {
      if(ibinx == 3) {
	lowBin = 3;
	highBin = 5;
      }
      if(ibinx == 4) {
	lowBin = 1;
	highBin = 5;
      }
      if(ibinx > 4) 
	continue;
    } //combineHigherJetMultBins
    
    double nWJets_ee	= GetEntries(h_wJets_FRee, lowBin, highBin);
    double nWJets_eeE	= GetTotalError(h_wJets_FRee, lowBin, highBin);
    double nQCD_ee	= GetEntries(h_qcd_FRee, lowBin, highBin);
    double nQCD_eeE	= GetTotalError(h_qcd_FRee, lowBin, highBin);
    
    double nWJets_mm	= GetEntries(h_wJets_FRmm, lowBin, highBin);
    double nWJets_mmE	= GetTotalError(h_wJets_FRmm, lowBin, highBin);
    double nQCD_mm	= GetEntries(h_qcd_FRmm, lowBin, highBin);
    double nQCD_mmE	= GetTotalError(h_qcd_FRmm, lowBin, highBin);

    double nWJets_em	= GetEntries(h_wJets_FRem, lowBin, highBin);
    double nWJets_emE	= GetTotalError(h_wJets_FRem, lowBin, highBin);
    double nQCD_em	= GetEntries(h_qcd_FRem, lowBin, highBin);
    double nQCD_emE	= GetTotalError(h_qcd_FRem, lowBin, highBin);

    double nWJets_all	= nWJets_ee + nWJets_mm + nWJets_em;
    double nWJets_allE	= sqrt(pow(nWJets_eeE,2) + pow(nWJets_mmE,2) + pow(nWJets_emE,2));
    double nQCD_all	= nQCD_ee + nQCD_mm + nQCD_em;
    double nQCD_allE	= sqrt(pow(nQCD_eeE,2) + pow(nQCD_mmE,2) + pow(nQCD_emE,2));


    
    double nFR_ee	=  nWJets_ee - nQCD_ee;
    double nFR_eeE	= sqrt(pow(nWJets_eeE,2) + pow(nQCD_eeE,2));

    double nFR_mm	=  nWJets_mm - nQCD_mm;
    double nFR_mmE	= sqrt(pow(nWJets_mmE,2) + pow(nQCD_mmE,2));

    double nFR_em	=  nWJets_em - nQCD_em;
    double nFR_emE	= sqrt(pow(nWJets_emE,2) + pow(nQCD_emE,2));

    double nFR_all	= nWJets_all - nQCD_all;
    double nFR_allE	= sqrt(pow(nWJets_allE,2) + pow(nQCD_allE,2));


    double nData_ee      = GetEntries(h_data_ee, lowBin, highBin); 
    double nData_mm      = GetEntries(h_data_mm, lowBin, highBin); 
    double nData_em      = GetEntries(h_data_em, lowBin, highBin); 
    double nData_all     = nData_ee + nData_mm + nData_em;

    double nData_eeE     = GetTotalError(h_data_ee, lowBin, highBin); 
    double nData_mmE     = GetTotalError(h_data_mm, lowBin, highBin); 
    double nData_emE     = GetTotalError(h_data_em, lowBin, highBin); 


    nWJets_ee_total = nWJets_ee_total + nWJets_ee;
    nWJets_eeE_total = sqrt(pow(nWJets_eeE_total,2) + pow(nWJets_eeE,2));
    nQCD_ee_total = nQCD_ee_total + nQCD_ee;
    nQCD_eeE_total = sqrt(pow(nQCD_eeE_total,2) + pow(nQCD_eeE,2));

    nWJets_mm_total = nWJets_mm_total + nWJets_mm;
    nWJets_mmE_total = sqrt(pow(nWJets_mmE_total,2) + pow(nWJets_mmE,2));
    nQCD_mm_total = nQCD_mm_total + nQCD_mm;
    nQCD_mmE_total = sqrt(pow(nQCD_mmE_total,2) + pow(nQCD_mmE,2));

    nWJets_em_total = nWJets_em_total + nWJets_em;
    nWJets_emE_total = sqrt(pow(nWJets_emE_total,2) + pow(nWJets_emE,2));
    nQCD_em_total = nQCD_em_total + nQCD_em;
    nQCD_emE_total = sqrt(pow(nQCD_emE_total,2) + pow(nQCD_emE,2));


    double SR_ee = nFR_ee/nData_ee;
    double SR_mm = nFR_mm/nData_mm;
    double SR_em = 0.5*(SR_ee + SR_mm);

    double SR_eeE = sqrt(pow(nFR_eeE/nData_ee, 2) + pow(nFR_ee*nData_eeE/pow(nData_ee,2), 2) );
    double SR_mmE = sqrt(pow(nFR_mmE/nData_mm, 2) + pow(nFR_mm*nData_mmE/pow(nData_mm,2), 2) );
    double SR_emE = 0.5*sqrt(pow(SR_mm*SR_eeE, 2) + pow(SR_ee*SR_mmE, 2));
    
    h_SR_ee->SetBinContent(ibinx, SR_ee);
    h_SR_ee->SetBinError(ibinx, SR_eeE);

    h_SR_mm->SetBinContent(ibinx, SR_mm);
    h_SR_mm->SetBinError(ibinx, SR_mmE);

    h_SR_em->SetBinContent(ibinx, SR_em);
    h_SR_em->SetBinError(ibinx, SR_emE);
    
        
    if(verbose) {
      
      if(combineHigherJetMultBins) {
	if(ibinx < 3)
	  cout << "\\textrm{" << ibinx-1 << " Jet Bin:}" << endL << endl;
	else if(ibinx == 3) 
	  cout << "\\textrm{$>=$ 2 Jet Bin:}" << endL << endl;
	else if(ibinx == 4) 
	  cout << "\\textrm{All Jets:} " << endL << endl;
	cout << "\\textrm{ } " << endL << endl;
      } else {
	  cout << "\\textrm{" << ibinx-1 << " Jet Bin:}" << endL << endl;
	  cout << "\\textrm{ } " << endL << endl;
      }

      cout << "\\begin{tabular}{l |  c  c  c  c} " << endl;
      cout << "\\hline" << endl;
      /*
	if(ibinx < 3)
	cout << ibinx-1 << " Jet Bin " << colSep << " " << colSep << " " << colSep << " " << colSep << endL << endl;
	else if(ibinx == 3)
	cout << "$>$ = 2 Jet Bin " << colSep << " " << colSep << " " << colSep << " " << colSep << endL << endl;
	else if(ibinx == 4)
	cout << "All Jets " << colSep << " " << colSep << " " << colSep << " " << colSep << endL << endl;
      */
      cout << Form("%20s", "Sample") << colSep << "ee" << colSep << "$\\mu\\mu$" << colSep << "e$\\mu$" << colSep << "all" << endL << " \\hline" << endl;
      cout << Form("%20s ", "WJets Prediction") << colSep 
	   << formatFloat(nWJets_ee, formatS)	<< pm	<< formatFloat(nWJets_eeE, formatS)	<< colSep
	   << formatFloat(nWJets_mm, formatS)	<< pm	<< formatFloat(nWJets_mmE, formatS)	<< colSep
	   << formatFloat(nWJets_em, formatS)	<< pm	<< formatFloat(nWJets_emE, formatS)	<< colSep
	   << formatFloat(nWJets_all, formatS)	<< pm	<< formatFloat(nWJets_allE, formatS)	<< endL << endl;
      cout << Form("%20s ", "QCD Prediction") << colSep 
	   << formatFloat(nQCD_ee, formatS) << pm << formatFloat(nQCD_eeE, formatS) << colSep
	   << formatFloat(nQCD_mm, formatS) << pm << formatFloat(nQCD_mmE, formatS) << colSep
	   << formatFloat(nQCD_em, formatS) << pm << formatFloat(nQCD_emE, formatS) << colSep
	   << formatFloat(nQCD_all, formatS) << pm << formatFloat(nQCD_allE, formatS) << endL << " \\hline" << endl;     
      cout << Form("%20s ", "FR Prediction") << colSep 
	   << formatFloat(nFR_ee,formatS) << pm << formatFloat(nFR_eeE,formatS) << colSep
	   << formatFloat(nFR_mm,formatS) << pm << formatFloat(nFR_mmE,formatS) << colSep
	   << formatFloat(nFR_em,formatS) << pm << formatFloat(nFR_emE,formatS) << colSep
	   << formatFloat(nFR_all,formatS) << pm << formatFloat(nFR_allE,formatS) << endL << endl;
      if(!haveAppliedZvetoCut) {
	if(useZmassRegionForSpillage)
	  cout << Form("%20s ", "Data Yield in Z mass") << colSep;
	else 
	  cout << Form("%20s ", "Data Yield") << colSep;
	cout << (int)nData_ee << colSep
	     << (int)nData_mm << colSep
	     << (int)nData_em << colSep 
	     << (int)nData_all << endL << endl;	
	cout << Form("%20s ", "Spillage Fraction") << colSep
	     << formatFloat(SR_ee, formatS) << pm << formatFloat(SR_eeE, formatS) << colSep 
	     << formatFloat(SR_mm, formatS) << pm << formatFloat(SR_mmE, formatS) << colSep
	     << formatFloat(SR_em, formatS) << pm << formatFloat(SR_mmE, formatS) << colSep
	     << " - " << colSep
	     << endL << " \\hline" << endl;
      }
      cout << "\\end{tabular}" << endl;
      cout << "\\vspace{2em} " << endL << endl;


    }


    
  }//loop
  
  f_FRwJets->Close();
  f_FRqcd->Close();
  
  v_SR.push_back(h_SR_ee);
  v_SR.push_back(h_SR_mm);
  v_SR.push_back(h_SR_em);

  return v_SR;
  
}

    
  



  

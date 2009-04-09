#include <iostream>
#include <vector>
#include <iomanip>

#include "TH1F.h"
#include "TH1F.h"
#include "TFile.h"
#include "TROOT.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"

using namespace std;


//void printR(TString fname, int nJets, TString hyp_type);
Double_t getError2(const TH1F *hist, const Int_t binMin, const Int_t binMax)
{

	//Double_t integral = 0.0;
	Double_t error2 = 0.0;
	for (Int_t i = binMin; i <= binMax; ++i)
	{
	  error2 += pow(hist->GetBinError(i), 2);
	}
	return error2;
}


pair<Double_t, Double_t> getR(TString fname, Int_t nJets, TString hyp_type) {

  TFile *f = TFile::Open(fname.Data(), "READ");
  
  //Get the DY dilMass histo for this jet bin and type
  TString histoname(Form("DY%s_hdilMass_%ij_%s", hyp_type.Data(), nJets,hyp_type.Data()));
  TH1F *h_mll = (TH1F*)f->Get(histoname.Data());
  if (h_mll==NULL) {
    cout << "Can't find " << histoname << ". Returning" << endl;
    return make_pair(0.,0.);
  }
  Int_t binLow = h_mll->FindBin(76.01);
  Int_t binHigh = h_mll->FindBin(105.09);
  Double_t nEvts_in = h_mll->Integral(binLow, binHigh);
  Double_t nEvts_out = h_mll->Integral(0, h_mll->GetNbinsX()+1) - nEvts_in;
  Double_t nEvts_inErr = sqrt(getError2(h_mll, binLow, binHigh));
  Double_t nEvts_outErr = sqrt(getError2(h_mll, 0, binLow) + getError2(h_mll, binHigh+1, h_mll->GetNbinsX()+1));
  float R = nEvts_out/(float)nEvts_in;
  //Error on R is R*R*(pow(nEvts_outErr/nEvts_out,2)+ pow(nEvts_inErr/nEvts_in,2))
  float R_Err2 = pow(nEvts_outErr/nEvts_in, 2) + pow(nEvts_out*nEvts_inErr,2)/pow(nEvts_in,4);
  
  
  cout << "R +/- err: " << R << "+/-" << sqrt(R_Err2) 
  << " (" << nEvts_in << "," << nEvts_out << ")" << endl;
  f->Close();
  return make_pair(R, R_Err2);

}
  
pair<Double_t, Double_t>  getkFromDY(TString fname, int nJets, TString hyp_type) {
  
  //getk in DY only (emu of course)
  TFile *f = TFile::Open(fname.Data(), "READ");
  TString histoname(Form("DYee_hdilMass_%ij_ee", nJets));
  TH1F *h_ee = (TH1F*)f->Get(histoname.Data());
  h_ee->TH1F::Sumw2();
  histoname = Form("DYmm_hdilMass_%ij_mm", nJets);
  TH1F *h_mm = (TH1F*)f->Get(histoname.Data());
  h_mm->TH1F::Sumw2();
  
  Int_t binLow = h_ee->FindBin(76.01);
  Int_t binHigh = h_ee->FindBin(105.09);
  Double_t nEvts_in_ee = h_ee->Integral(binLow, binHigh);
  Double_t nEvts_in_eeErr2 = getError2(h_ee, binLow, binHigh);

  binLow = h_mm->FindBin(76.01);
  binHigh = h_mm->FindBin(105.09);
  Double_t nEvts_in_mm = h_mm->Integral(binLow, binHigh);
  Double_t nEvts_in_mmErr2 = getError2(h_mm, binLow, binHigh);
  
  Double_t k, kErr2;
  
  if(hyp_type == "mm") {
    k = 0.5*sqrt(nEvts_in_mm/nEvts_in_ee);
    kErr2 = (1./16.)*(nEvts_in_mmErr2/nEvts_in_ee/nEvts_in_mm + nEvts_in_eeErr2*nEvts_in_mm/pow(nEvts_in_ee,3));
  } else if(hyp_type == "ee" ) {
    k = 0.5*sqrt(nEvts_in_ee/nEvts_in_mm);
    kErr2 = (1./16.)*(nEvts_in_eeErr2/nEvts_in_mm/nEvts_in_ee + nEvts_in_mmErr2*nEvts_in_ee/pow(nEvts_in_mm,3));
  }
  //cout << "k " << k << " +/- " << kErr << endl;

  //Double_t eff_movere = nEvts_in_mm/nEvts_in_ee;
  //Double_t eff_movere_Err = sqrt(nEvts_in_mmErr2/pow(nEvts_in_ee,2) + pow(nEvts_in_mm,2)*nEvts_in_eeErr2/pow(nEvts_in_ee,4));

  //cout << "efficiency, " << nJets << "jets: " << eff_movere << " +/- " << eff_movere_Err << endl;
  
  f->Close();
  return make_pair(k, kErr2);
  
}
  

pair<Double_t, Double_t> getkFromAll(TString fname, int nJets, TString hyp_type) {
  
  
  TFile *f = TFile::Open(fname.Data(), "READ");
  vector<TH1F*> vhistos_ee;
  vector<TH1F*> vhistos_mm;
  
  TList *l_keys = f->GetListOfKeys();
  for(int i = 0; i< l_keys->GetSize(); i++) {
    TString histoName(l_keys->At(i)->GetName());
    if(histoName.Contains("WQQ")) continue;
    if(histoName.Contains("wz")) continue;
    if(histoName.EndsWith("ee") && histoName.Contains(Form("hdilMass_%ij_",nJets) ) ) 
      vhistos_ee.push_back((TH1F*)f->Get(histoName.Data()));
    if(histoName.EndsWith("mm") && histoName.Contains(Form("hdilMass_%ij_",nJets) ) ) 
      vhistos_mm.push_back((TH1F*)f->Get(histoName.Data()));
  }


  //get the numbers for ee
  Double_t nEvts_in_ee = 0.0;
  Double_t nEvts_in_eeErr2 = 0.0;
  for(unsigned int i = 0 ; i < vhistos_ee.size(); i++) {
    Int_t binLow  = vhistos_ee.at(i)->FindBin(76.01);
    Int_t binHigh = vhistos_ee.at(i)->FindBin(105.09);
    nEvts_in_ee += vhistos_ee.at(i)->Integral(binLow, binHigh);
    nEvts_in_eeErr2 += getError2(vhistos_ee.at(i), binLow, binHigh);
  }
  
  Double_t nEvts_in_mm = 0.0;
  Double_t nEvts_in_mmErr2 = 0.0;
  for(unsigned int i = 0 ; i < vhistos_mm.size(); i++) {
    Int_t binLow  = vhistos_mm.at(i)->FindBin(76.01);
    Int_t binHigh = vhistos_mm.at(i)->FindBin(105.09);
    nEvts_in_mm += vhistos_mm.at(i)->Integral(binLow, binHigh);
    nEvts_in_mmErr2 += getError2(vhistos_mm.at(i), binLow, binHigh);
  }
  
  Double_t k, kErr2;
  if(hyp_type == "mm") {
    k = 0.5*sqrt(nEvts_in_mm/nEvts_in_ee);
    kErr2 = (1./16.)*(nEvts_in_mmErr2/nEvts_in_ee/nEvts_in_mm + nEvts_in_eeErr2*nEvts_in_mm/pow(nEvts_in_ee,3));
  } else if(hyp_type == "ee") {
    k = 0.5*sqrt(nEvts_in_ee/nEvts_in_mm);
    kErr2 = (1./16.)*(nEvts_in_eeErr2/nEvts_in_mm/nEvts_in_ee + nEvts_in_mmErr2*nEvts_in_ee/pow(nEvts_in_mm,3));
  }
    
    
    
  //cout << "k " << k << " +/- " << kErr << endl;
  
  //Double_t eff_movere = nEvts_in_mm/nEvts_in_ee;
  //Double_t eff_movere_Err2 = nEvts_in_mmErr2/pow(nEvts_in_ee,2) + pow(nEvts_in_mm,2)*nEvts_in_eeErr2/pow(nEvts_in_ee,4);
  
  //cout << "efficiency, " << nJets << "jets: " << eff_movere << " +/- " << eff_movere_Err << endl;
  f->Close();
  return make_pair(k, kErr2);
}
    
  
pair<Double_t,Double_t> getTrueNOut(TString fname, int nJets, TString hyptype ) {

  
  TFile *f = TFile::Open(fname.Data(), "READ");
  
  //Get the DY dilMass histo for this jet bin and type
  TString histoname(Form("DY%s_hdilMass_%ij_%s", hyptype.Data(), nJets,hyptype.Data()));
  TH1F *h_dy_mll = (TH1F*)f->Get(histoname.Data());
  histoname = Form("zz_hdilMass_%ij_%s", nJets,hyptype.Data());
  TH1F *h_zz_mll = (TH1F*)f->Get(histoname.Data());
  if (h_dy_mll==NULL) {
    cout << "Can't find " << histoname << ". Returning" << endl;
    return make_pair(0.,0.);
  }

  if (h_zz_mll==NULL) {
    cout << "Can't find " << histoname << ". Returning" << endl;
    return make_pair(0.,0.);
  }
    
  Int_t binLow  = h_dy_mll->FindBin(76.01);
  Int_t binHigh = h_dy_mll->FindBin(105.09);
  Double_t nOut = h_dy_mll->Integral(0, binLow-1) + h_dy_mll->Integral(binHigh+1, h_dy_mll->GetNbinsX()+1);
  Double_t nOutErr2 = getError2(h_dy_mll, 0, binLow-1) + getError2(h_dy_mll, binHigh+1, h_dy_mll->GetNbinsX()+1);

  binLow  = h_zz_mll->FindBin(76.01);
  binHigh = h_zz_mll->FindBin(105.09);
  nOut = nOut + h_zz_mll->Integral(0, binLow-1) + h_zz_mll->Integral(binHigh+1, h_zz_mll->GetNbinsX()+1);
  nOutErr2 = nOutErr2 + getError2(h_zz_mll, 0, binLow-1) + getError2(h_zz_mll, binHigh+1, h_zz_mll->GetNbinsX()+1);
  
  delete h_dy_mll;
  delete h_zz_mll;
  
  f->Close();
  return make_pair(nOut, nOutErr2);

}

pair<Double_t, Double_t> getnEMuIn(TString fname, int nJets) {
 
  TFile *f = TFile::Open(fname.Data(), "READ");
  vector<TH1F*> vhistos_em;
  
  TList *l_keys = f->GetListOfKeys();
  for(int i = 0; i< l_keys->GetSize(); i++) {
    TString histoName(l_keys->At(i)->GetName());
    if(histoName.Contains("WQQ")) continue;
    if(histoName.Contains("wz")) continue;
    if(histoName.EndsWith("em") && histoName.Contains(Form("hdilMass_%ij_",nJets) ) ) 
      vhistos_em.push_back((TH1F*)f->Get(histoName.Data()));
    }

  
  Double_t nEMuIn = 0.;
  Double_t nEMuInErr2 = 0.;
  for(unsigned int i = 0; i < vhistos_em.size(); i++) {
    Int_t binLow  = vhistos_em.at(i)->FindBin(76.01);
    Int_t binHigh = vhistos_em.at(i)->FindBin(105.09);
    nEMuIn += vhistos_em.at(i)->Integral(binLow, binHigh);
    nEMuInErr2 = getError2(vhistos_em.at(i), binLow, binHigh);
  }
  
  return make_pair(nEMuIn, nEMuInErr2 );  
  f->Close();
}


pair<Double_t, Double_t> getnLLIn(TString fname, int nJets, TString hyptype) {
 
  TFile *f = TFile::Open(fname.Data(), "READ");
  vector<TH1F*> vhistos;
  
  TList *l_keys = f->GetListOfKeys();
  for(int i = 0; i< l_keys->GetSize(); i++) {
    TString histoName(l_keys->At(i)->GetName());
    if(histoName.Contains("WQQ")) continue;
    if(histoName.Contains("wz")) continue;
    if(histoName.Contains(Form("hdilMass_%ij_%s",nJets, hyptype.Data()) ) )
      vhistos.push_back((TH1F*)f->Get(histoName.Data()));
    }

  
  Double_t nllIn = 0.;
  Double_t nllInErr2 = 0.;
  for(unsigned int i = 0; i < vhistos.size(); i++) {
    Int_t binLow  = vhistos.at(i)->FindBin(76.01);
    Int_t binHigh = vhistos.at(i)->FindBin(105.09);
    nllIn += vhistos.at(i)->Integral(binLow, binHigh);
    nllInErr2 = getError2(vhistos.at(i), binLow, binHigh);
  }
  
  return make_pair(nllIn, nllInErr2 );  
  f->Close();
  vhistos.clear();
}



void getPredicted() {
  
  //TString f_NoMETCut("myHist_isoReg_looseDil08_OS_noDupWt_isoDil08_hltMu9E15.root");
  //TString f_wMETCut("myHist_isoReg_looseDil08_OS_noDupWt_isoDil08_preMet08_hltMu9E15.root");
  
  TString f_NoMETCut("myHist_1106944__looseDil08_OS_noDupWt_isoDil08_hltMu9E15.root");
  TString f_wMETCut("myHist_1369088__looseDil08_OS_noDupWt_isoDil08_preMet08_hltMu9E15.root");

  //Declare histos for predicted numbers
  TH1F *h_Nout_ee = new TH1F("h_Nout_ee", "", 3, 0,3);
  TH1F *h_Nout_eeTrue = new TH1F("h_Nout_ee", "", 3, 0,3);
  TH1F *h_Nout_eeTrueH = new TH1F("h_Nout_eeH", "", 3, 0,3);
  TH1F *h_Nout_eeTrueL = new TH1F("h_Nout_eeL", "", 3, 0,3);

  TH1F *h_Nout_mm = new TH1F("h_Nout_mm", "", 3, 0,3);
  TH1F *h_Nout_mmTrue = new TH1F("h_Nout_mm", "", 3, 0,3);
  TH1F *h_Nout_mmTrueH = new TH1F("h_Nout_mmH", "", 3, 0,3);
  TH1F *h_Nout_mmTrueL = new TH1F("h_Nout_mmL", "", 3, 0,3);

  
  TAxis *xAxis_ee = h_Nout_ee->GetXaxis();
  xAxis_ee->SetTitle("nJets");
  xAxis_ee->SetBinLabel(1,"0");
  xAxis_ee->SetBinLabel(2,"1");
  xAxis_ee->SetBinLabel(3,">=2");
  TAxis *xAxis_eeTrue = h_Nout_eeTrue->GetXaxis();
  xAxis_eeTrue->SetTitle("nJets");
  xAxis_eeTrue->SetBinLabel(1,"0");
  xAxis_eeTrue->SetBinLabel(2,"1");
  xAxis_eeTrue->SetBinLabel(3,">=2");
  TAxis *xAxis_eeTrueH = h_Nout_eeTrueH->GetXaxis();
  xAxis_eeTrueH->SetTitle("nJets");
  xAxis_eeTrueH->SetBinLabel(1,"0");
  xAxis_eeTrueH->SetBinLabel(2,"1");
  xAxis_eeTrueH->SetBinLabel(3,">=2");
  TAxis *xAxis_eeTrueL = h_Nout_eeTrueL->GetXaxis();
  xAxis_eeTrueL->SetTitle("nJets");
  xAxis_eeTrueL->SetBinLabel(1,"0");
  xAxis_eeTrueL->SetBinLabel(2,"1");
  xAxis_eeTrueL->SetBinLabel(3,">=2");

  TAxis *xAxis_mm = h_Nout_mm->GetXaxis();
  xAxis_mm->SetTitle("nJets");
  xAxis_mm->SetBinLabel(1,"0");
  xAxis_mm->SetBinLabel(2,"1");
  xAxis_mm->SetBinLabel(3,">=2");
  TAxis *xAxis_mmTrue = h_Nout_mmTrue->GetXaxis();
  xAxis_mmTrue->SetTitle("nJets");
  xAxis_mmTrue->SetBinLabel(1,"0");
  xAxis_mmTrue->SetBinLabel(2,"1");
  xAxis_mmTrue->SetBinLabel(3,">=2");
  TAxis *xAxis_mmTrueH = h_Nout_mmTrueH->GetXaxis();
  xAxis_mmTrueH->SetTitle("nJets");
  xAxis_mmTrueH->SetBinLabel(1,"0");
  xAxis_mmTrueH->SetBinLabel(2,"1");
  xAxis_mmTrueH->SetBinLabel(3,">=2");
  TAxis *xAxis_mmTrueL = h_Nout_mmTrueL->GetXaxis();
  xAxis_mmTrueL->SetTitle("nJets");
  xAxis_mmTrueL->SetBinLabel(1,"0");
  xAxis_mmTrueL->SetBinLabel(2,"1");
  xAxis_mmTrueL->SetBinLabel(3,">=2");
    
  
  vector<pair<Double_t, Double_t> > v_k_ee;
  vector<pair<Double_t, Double_t> > v_k_mm;
  
  //ee case
  cout << "EE" << endl;
  cout << "|  *nJets*  |  *NIn_ll*  |  *NIn_emu*  |  *R*  |  *Nout(Predicted)*  |  *Nout(True)*  |" << endl;
  
  for(int nJet = 0; nJet < 3; nJet++) {
    //get Rout/Rin
    pair<Double_t, Double_t> p_rOutIn_ee = getR(f_wMETCut, nJet, "ee");
    
    //get k 
    pair<Double_t, Double_t> p_k = getkFromAll(f_NoMETCut, nJet, "ee");
    v_k_ee.push_back(p_k);
    //get the true number of events in the out region in DY
    pair<Double_t, Double_t> p_truNOut_ee = getTrueNOut(f_wMETCut, nJet, "ee");
    
    //get the number of events in the mass range, emu
    pair<Double_t, Double_t> p_nIn_emu = getnEMuIn(f_wMETCut, nJet);

    //get the number of events in the dilepton mass range
    pair<Double_t, Double_t> p_nllIn_ee = getnLLIn(f_wMETCut, nJet, "ee");
    
    Double_t predicted_ee = (p_nllIn_ee.first - p_k.first*p_nIn_emu.first)*p_rOutIn_ee.first;
    
    Double_t dpredict_dR      = p_nllIn_ee.first-p_k.first*p_nIn_emu.first;
    Double_t dpredict_dk      = -p_rOutIn_ee.first*p_nIn_emu.first;
    Double_t dpredict_dnllin  = p_rOutIn_ee.first;
    Double_t dpredict_dnInemu = -p_rOutIn_ee.first*p_k.first;

    Double_t rOutInErr2_ee  = p_rOutIn_ee.second;
    Double_t kErr2_ee       = p_k.second;
    Double_t nllInErr2_ee   = p_nllIn_ee.second;
    Double_t nemuInErr2_ee  = p_nIn_emu.second;
    

    Double_t predictedErr2_ee = pow(dpredict_dR,2)*rOutInErr2_ee + pow(dpredict_dk,2)*kErr2_ee +
      pow(dpredict_dnllin,2)*nllInErr2_ee + pow(dpredict_dnInemu,2)*nemuInErr2_ee;
      
    
    h_Nout_ee->SetBinContent(nJet+1, predicted_ee);
    h_Nout_ee->SetBinError(nJet+1, sqrt(predictedErr2_ee));
    h_Nout_eeTrue->SetBinContent(nJet+1, p_truNOut_ee.first);
    h_Nout_eeTrueH->SetBinContent(nJet+1, p_truNOut_ee.first + sqrt(p_truNOut_ee.second));
    h_Nout_eeTrueL->SetBinContent(nJet+1, p_truNOut_ee.first - sqrt(p_truNOut_ee.second));
    
    
    cout <<"|  " <<  nJet << "  |  "
      	 << p_nllIn_ee.first << " +/- " << sqrt(p_nllIn_ee.second) << "  |  " 
	 << p_nIn_emu.first << " +/- " << sqrt(p_nIn_emu.second) << "  |  "
	 << p_rOutIn_ee.first << " +/- " << sqrt(p_rOutIn_ee.second) << "  |  "
	 << predicted_ee << " +/- " << sqrt(predictedErr2_ee) << "  |  "
	 << p_truNOut_ee.first << " +/- " << sqrt(p_truNOut_ee.second) << "  |  " << endl;

    
    Double_t nDYErr2 = p_nllIn_ee.first + pow(p_k.first,2)*p_nIn_emu.first + pow(p_nIn_emu.first,2)*kErr2_ee;
    //Double_t dlEvans_err = pow(predicted_ee,2)*(rOutInErr2_ee/pow(p_rOutIn_ee.first,2) + nDYErr2/pow(p_nllIn_ee.first,2));
    
    //cout << sqrt(dlEvans_err) << endl;
    
  }
  
  cout << endl << "| *nJets*  |  *k*  |" << endl;
  for(unsigned int i = 0 ; i < v_k_ee.size(); i++) 
    cout << "|  " << i << "  |  " 
	 << v_k_ee.at(i).first << " +/- " << sqrt(v_k_ee.at(i).second) << "  |" << endl;
  
  
  cout << endl << endl << "MM" << endl;
  //cout << "|  *nJets*  |  *R*  |  *Nout(Predicted)*  |  *Nout(True)*  |" << endl;
  cout << "|  *nJets*  |  *NIn_ll*  |  *NIn_emu*  |  *R*  |  *Nout(Predicted)*  |  *Nout(True)*  |" << endl;
  
  //mm
  for(int nJet = 0; nJet < 3; nJet++) {
    //get Rout/Rin
    pair<Double_t, Double_t> p_rOutIn_mm = getR(f_wMETCut, nJet, "mm");
    
    //get k 
    pair<Double_t, Double_t> p_k = getkFromAll(f_NoMETCut, nJet, "mm");
    v_k_mm.push_back(p_k);
    
    //get the true number of events in the out region in DY
    pair<Double_t, Double_t> p_truNOut_mm = getTrueNOut(f_wMETCut, nJet, "mm");

    //get the number of events in the mass range, emu
    pair<Double_t, Double_t> p_nIn_emu = getnEMuIn(f_wMETCut, nJet);

    //get the number of events in the dilepton mass range
    pair<Double_t, Double_t> p_nllIn_mm = getnLLIn(f_wMETCut, nJet, "mm");
    
    Double_t predicted_mm = (p_nllIn_mm.first - p_k.first*p_nIn_emu.first)*p_rOutIn_mm.first;

 
    Double_t dpredict_dR      = p_nllIn_mm.first-p_k.first*p_nIn_emu.first;
    Double_t dpredict_dk      = -p_rOutIn_mm.first*p_nIn_emu.first;
    Double_t dpredict_dnllin  = p_rOutIn_mm.first;
    Double_t dpredict_dnInemu = -p_rOutIn_mm.first*p_k.first;

    Double_t rOutInErr2_mm  = p_rOutIn_mm.second;
    Double_t kErr2_mm       = p_k.second;
    Double_t nllInErr2_mm   = p_nllIn_mm.second;
    Double_t nemuInErr2_mm  = p_nIn_emu.second;
    

    Double_t predictedErr2_mm = pow(dpredict_dR,2)*rOutInErr2_mm + pow(dpredict_dk,2)*kErr2_mm +
      pow(dpredict_dnllin,2)*nllInErr2_mm + pow(dpredict_dnInemu,2)*nemuInErr2_mm;
      
      
    //mm
    cout << "|  " << nJet << "  |  " 
	 << p_nllIn_mm.first << " +/- " << sqrt(p_nllIn_mm.second) << "  |  "
	 << p_nIn_emu.first << " +/- " << sqrt(p_nIn_emu.second) << "  |  "
	 << p_rOutIn_mm.first << " +/- " << sqrt(p_rOutIn_mm.second) << "  |  "	 
	 << predicted_mm << " +/- " << sqrt(predictedErr2_mm) << "  |  "
	 << p_truNOut_mm.first << " +/- " << sqrt(p_truNOut_mm.second) << "  |" << endl;
    
    Double_t nDYErr2 = p_nllIn_mm.first + pow(p_k.first,2)*p_nIn_emu.first + pow(p_nIn_emu.first,2)*kErr2_mm;
    //Double_t dlEvans_err = pow(predicted_mm,2)*(rOutInErr2_mm/pow(p_rOutIn_mm.first,2) + nDYErr2/pow(p_nllIn_mm.first,2));
    
    //cout << sqrt(dlEvans_err) << endl;
    
    h_Nout_mm->SetBinContent(nJet+1, predicted_mm);
    h_Nout_mm->SetBinError(nJet+1, sqrt(predictedErr2_mm));
    h_Nout_mmTrue->SetBinContent(nJet+1, p_truNOut_mm.first);
    h_Nout_mmTrueH->SetBinContent(nJet+1, p_truNOut_mm.first + sqrt(p_truNOut_mm.second));
    h_Nout_mmTrueL->SetBinContent(nJet+1, p_truNOut_mm.first - sqrt(p_truNOut_mm.second));

    
  }

  gStyle->SetOptStat(0);
  cout << endl << "| *nJets*  |  *k*  |" << endl;
  for(unsigned int i = 0 ; i < v_k_mm.size(); i++) 
    cout << "|  " << i << "  |  " 
	 << v_k_mm.at(i).first << " +/- " << sqrt(v_k_mm.at(i).second) << "  |" << endl;
  
  TCanvas *cee = new TCanvas();
  cee->cd();
  
  h_Nout_eeTrueH->GetYaxis()->SetRangeUser(0,10);
  h_Nout_eeTrueH->SetFillColor(kYellow);
  h_Nout_eeTrueH->SetLineColor(10);
  h_Nout_eeTrueH->Draw("hist");
  
  h_Nout_eeTrueL->GetYaxis()->SetRangeUser(0,10);
  h_Nout_eeTrueL->SetFillColor(10);
  h_Nout_eeTrueL->SetLineColor(10);
  h_Nout_eeTrueL->SetFillStyle(1001);
  h_Nout_eeTrueL->Draw("histsame");  

  h_Nout_eeTrue->GetYaxis()->SetRangeUser(0,10);
  h_Nout_eeTrue->SetMinimum(0.9*h_Nout_eeTrueL->GetMinimum());
  h_Nout_eeTrue->SetMaximum(1.1*h_Nout_eeTrueH->GetMaximum());
  h_Nout_eeTrue->Draw("same");  

  h_Nout_ee->GetYaxis()->SetRangeUser(0,10);
  h_Nout_ee->SetMinimum(0.9*h_Nout_eeTrueL->GetMinimum());
  h_Nout_ee->SetMaximum(1.1*h_Nout_eeTrueH->GetMaximum());
  h_Nout_ee->SetMarkerStyle(20);
  //h_Nout_ee->SetMarkerSize(2);
  h_Nout_ee->Draw("samee");  
  
  cee->RedrawAxis();

  TCanvas *cmm = new TCanvas();
  cmm->cd();
  
  h_Nout_mmTrueH->GetYaxis()->SetRangeUser(0,10);
  h_Nout_mmTrueH->SetFillColor(kYellow);
  h_Nout_mmTrueH->SetLineColor(10);
  h_Nout_mmTrueH->Draw("hist");
  
  h_Nout_mmTrueL->GetYaxis()->SetRangeUser(0,10);
  h_Nout_mmTrueL->SetFillColor(10);
  h_Nout_mmTrueL->SetLineColor(10);
  h_Nout_mmTrueL->SetFillStyle(1001);
  h_Nout_mmTrueL->Draw("histsame");  

  h_Nout_mmTrue->GetYaxis()->SetRangeUser(0,10);
  h_Nout_mmTrue->SetMinimum(0.9*h_Nout_mmTrueL->GetMinimum());
  h_Nout_mmTrue->SetMaximum(1.1*h_Nout_mmTrueH->GetMaximum());
  h_Nout_mmTrue->Draw("same");  

  h_Nout_mm->GetYaxis()->SetRangeUser(0,10);
  h_Nout_mm->SetMinimum(0.9*h_Nout_mmTrueL->GetMinimum());
  h_Nout_mm->SetMaximum(1.1*h_Nout_mmTrueH->GetMaximum());
  h_Nout_mm->SetMarkerStyle(20);
  //h_Nout_mm->SetMarkerSize(2);
  h_Nout_mm->Draw("samee");  
  
  cmm->RedrawAxis();

  

}
    

    
    
  

  

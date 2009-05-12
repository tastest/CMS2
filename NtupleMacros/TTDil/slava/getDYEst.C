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
  TString histoname_in(Form("DY%s_hnJetinZwindow_%s", hyp_type.Data(), hyp_type.Data()));
  TString histoname_out(Form("DY%s_hnJetoutZwindow_%s", hyp_type.Data(), hyp_type.Data()));
  TH1F *h_nJets_in = (TH1F*)f->Get(histoname_in.Data());
  TH1F *h_nJets_out = (TH1F*)f->Get(histoname_out.Data());
  
  
  Double_t nEvts_in  = 0;
  Double_t nEvts_out = 0;
  Double_t nEvts_inErr  = 0;
  Double_t nEvts_outErr = 0;

  if(nJets <= 1) {
    nEvts_in  = h_nJets_in->GetBinContent(nJets + 1);
    nEvts_out = h_nJets_out->GetBinContent(nJets + 1);
    nEvts_inErr  = sqrt(getError2(h_nJets_in, nJets + 1, nJets + 1));
    nEvts_outErr = sqrt(getError2(h_nJets_out, nJets + 1, nJets + 1));
  } else {
    Int_t binHigh = h_nJets_in->GetNbinsX() + 1;
    for(Int_t i = 3; i <= binHigh; i++) {
      nEvts_in = nEvts_in + h_nJets_in->GetBinContent(i);
      nEvts_out = nEvts_out + h_nJets_out->GetBinContent(i);
    }
    nEvts_inErr = sqrt(getError2(h_nJets_in, 3, binHigh));
    nEvts_outErr = sqrt(getError2(h_nJets_out, 3, binHigh));
  }
  
  Double_t R = nEvts_out/nEvts_in;
  Double_t R_Err2 = pow(nEvts_outErr/nEvts_in, 2) + pow(nEvts_out*nEvts_inErr,2)/pow(nEvts_in,4);
  
  f->Close();
  return make_pair(R, R_Err2);

}
  
pair<Double_t, Double_t>  getkFromDY(TString fname, int nJets, TString hyp_type) {
  
  //getk in DY only (emu of course)
   TFile *f = TFile::Open(fname.Data(), "READ");
   TString histoname("DYee_hnJetinZwindow_ee");
   TH1F *h_ee = (TH1F*)f->Get(histoname.Data());
   h_ee->TH1F::Sumw2();
   histoname = "DYmm_hnJetinZwindow_mm";
   TH1F *h_mm = (TH1F*)f->Get(histoname.Data());
   h_mm->TH1F::Sumw2();
   
   Double_t nEvts_in_ee     = 0;
   Double_t nEvts_inErr2_ee  = 0;
   Double_t nEvts_in_mm     = 0;
   Double_t nEvts_inErr2_mm  = 0;
   
   if(nJets <= 1) {
     nEvts_in_ee  = h_ee->GetBinContent(nJets+1);
     nEvts_in_mm  = h_mm->GetBinContent(nJets+1);
     nEvts_inErr2_ee = getError2(h_ee, nJets+1, nJets+1);
     nEvts_inErr2_mm = getError2(h_mm, nJets+1, nJets+1);
   } else {
     Int_t binHigh = h_ee->GetNbinsX() + 1;
     for(Int_t i = 3; i <= binHigh; i++) {
       nEvts_in_ee = nEvts_in_ee + h_ee->GetBinContent(i);
       nEvts_in_mm = nEvts_in_mm + h_ee->GetBinContent(i);
     }
     nEvts_inErr2_ee = getError2(h_ee, 3, binHigh);
     nEvts_inErr2_mm = getError2(h_mm, 3, binHigh);
   }
			   
   Double_t k     = 0.;
   Double_t kErr2 = 0.;
   if(hyp_type == "mm") {
      k = 0.5*sqrt(nEvts_in_mm/nEvts_in_ee);
      kErr2 = (1./16.)*(nEvts_inErr2_mm/nEvts_in_ee/nEvts_in_mm + nEvts_inErr2_ee*nEvts_in_mm/pow(nEvts_in_ee,3));
   } else if(hyp_type == "ee" ) {
    k = 0.5*sqrt(nEvts_in_ee/nEvts_in_mm);
    kErr2 = (1./16.)*(nEvts_inErr2_ee/nEvts_in_mm/nEvts_in_ee + nEvts_inErr2_mm*nEvts_in_ee/pow(nEvts_in_mm,3));
   }

   f->Close();
   return make_pair(k,kErr2);
  
}
  

pair<Double_t, Double_t> getkFromAll(TString fname, int nJets, TString hyp_type) {
  
  
  TFile *f = TFile::Open(fname.Data(), "READ");
  vector<TH1F*> vhistos_ee;
  vector<TH1F*> vhistos_mm;
  
  TList *l_keys = f->GetListOfKeys();
  for(int i = 0; i< l_keys->GetSize(); i++) {
    TString histoName(l_keys->At(i)->GetName());
    if(histoName.EndsWith("ee") && histoName.Contains("hnJetinZwindow"))
      vhistos_ee.push_back((TH1F*)f->Get(histoName.Data()));
    if(histoName.EndsWith("mm") && histoName.Contains("hnJetinZwindow"))
      vhistos_mm.push_back((TH1F*)f->Get(histoName.Data()));
  }
  
 //get the numbers for ee
  Double_t nEvts_in_ee = 0.0;
  Double_t nEvts_inErr2_ee = 0.0;
  for(unsigned int j = 0 ; j < vhistos_ee.size(); j++) {
    if(nJets <= 1) {
      nEvts_in_ee     += vhistos_ee.at(j)->GetBinContent(nJets+1);
      nEvts_inErr2_ee += getError2(vhistos_ee.at(j), nJets+1, nJets+1);
    } else {
      Int_t binHigh = vhistos_ee.at(j)->GetNbinsX() + 1;
      for(Int_t i = 3; i <= binHigh; i++) 
	nEvts_in_ee +=vhistos_ee.at(j)->GetBinContent(i);
      nEvts_inErr2_ee = getError2(vhistos_ee.at(j), 3, binHigh);
    }
  }
  
  Double_t nEvts_in_mm = 0.0;
  Double_t nEvts_inErr2_mm = 0.0;
  for(unsigned int j = 0 ; j < vhistos_mm.size(); j++) {
    if(nJets <=1) {
      nEvts_in_mm += vhistos_mm.at(j)->GetBinContent(nJets+1);
      nEvts_inErr2_mm += getError2(vhistos_mm.at(j), nJets+1, nJets+1);
    } else {
      Int_t binHigh = vhistos_mm.at(j)->GetNbinsX() + 1;
      for(Int_t i = 3; i <= binHigh; i++) 
	nEvts_in_mm +=   vhistos_mm.at(j)->GetBinContent(i);
      nEvts_inErr2_mm += getError2(vhistos_mm.at(j), 3, binHigh);
    }
  }
  
  Double_t k     = 0.;
  Double_t kErr2 = 0.; 
  if(hyp_type == "mm") {
    k = 0.5*sqrt(nEvts_in_mm/nEvts_in_ee);
    kErr2 = (1./16.)*(nEvts_inErr2_mm/nEvts_in_ee/nEvts_in_mm + nEvts_inErr2_ee*nEvts_in_mm/pow(nEvts_in_ee,3));
  } else if(hyp_type == "ee") {
    k = 0.5*sqrt(nEvts_in_ee/nEvts_in_mm);
    kErr2 = (1./16.)*(nEvts_inErr2_ee/nEvts_in_mm/nEvts_in_ee + nEvts_inErr2_mm*nEvts_in_ee/pow(nEvts_in_mm,3));
  }

  f->Close();
  return make_pair(k, kErr2);

}
    
  
pair<Double_t,Double_t> getTrueNOut(TString fname, int nJets, TString hyptype ) {

  
  TFile *f = TFile::Open(fname.Data(), "READ");
  
  //Get the DY dilMass histo for this jet bin and type
  TString histoname(Form("DY%s_hnJetoutZwindow_%s", hyptype.Data(),hyptype.Data()));
  TH1F *h_nJets_dy = (TH1F*)f->Get(histoname.Data());
  histoname = Form("zz_hnJetoutZwindow_%s", hyptype.Data());
  TH1F *h_nJets_zz = (TH1F*)f->Get(histoname.Data());
  if (h_nJets_dy==NULL) {
    cout << "Can't find " << histoname << ". Returning" << endl;
    return make_pair(0.,0.);
  }

  if (h_nJets_zz==NULL) {
    cout << "Can't find " << histoname << ". Returning" << endl;
    return make_pair(0.,0.);
  }
    
  Double_t nOut     = 0.;
  Double_t nOutErr2 = 0.;
  
  if(nJets <=1) {
    nOut = h_nJets_dy->GetBinContent(nJets+1) + h_nJets_zz->GetBinContent(nJets+1);
    nOutErr2 = getError2(h_nJets_dy, nJets+1, nJets+1) + getError2(h_nJets_zz, nJets+1, nJets+1);
  } else {
    Int_t binHigh = h_nJets_dy->GetNbinsX()+1;
    for(Int_t i = 3; i <= binHigh; i++) 
      nOut     += h_nJets_dy->GetBinContent(i) + h_nJets_zz->GetBinContent(i);
    nOutErr2 += getError2(h_nJets_dy, 3, binHigh) + getError2(h_nJets_zz, 3, binHigh);
  }

  delete h_nJets_dy;
  delete h_nJets_zz;

  f->Close();
  return make_pair(nOut, nOutErr2);


}

pair<Double_t, Double_t> getnEMuIn(TString fname, int nJets) {
 
  TFile *f = TFile::Open(fname.Data(), "READ");
  vector<TH1F*> vhistos_em;
  
  TList *l_keys = f->GetListOfKeys();
  for(int i = 0; i< l_keys->GetSize(); i++) {
    TString histoName(l_keys->At(i)->GetName());
    if(histoName.EndsWith("em") && histoName.Contains("hnJetinZwindow"))
      vhistos_em.push_back((TH1F*)f->Get(histoName.Data()));
    }

  
  Double_t nEMuIn = 0.;
  Double_t nEMuInErr2 = 0.;
  for(unsigned int j = 0; j < vhistos_em.size(); j++) {
    if(nJets <= 1) {
      nEMuIn      += vhistos_em.at(j)->GetBinContent(nJets+1);
      nEMuInErr2  += getError2(vhistos_em.at(j), nJets+1, nJets+1);
    } else {
      Int_t binHigh = vhistos_em.at(j)->GetNbinsX() + 1;
      for(Int_t i = 3; i <= binHigh; i++) 
	nEMuIn   += vhistos_em.at(j)->GetBinContent(i);
      nEMuInErr2 += getError2(vhistos_em.at(j), 3, binHigh);
    }
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
    if(histoName.Contains(Form("hnJetinZwindow_%s", hyptype.Data()) ) )
      vhistos.push_back((TH1F*)f->Get(histoName.Data()));
    }

  
  Double_t nllIn = 0.;
  Double_t nllInErr2 = 0.;
  for(unsigned int j = 0; j < vhistos.size(); j++) {
    if(nJets <= 1) {
      nllIn      += vhistos.at(j)->GetBinContent(nJets+1);
      nllInErr2  += getError2(vhistos.at(j), nJets+1, nJets+1);
    } else {
      Int_t binHigh = vhistos.at(j)->GetNbinsX() + 1;
      for(Int_t i = 3; i <= binHigh; i++) 
	nllIn   += vhistos.at(j)->GetBinContent(i);
      nllInErr2 += getError2(vhistos.at(j), 3, binHigh);
    }
  }
  
  return make_pair(nllIn, nllInErr2 );  
  f->Close();
  vhistos.clear();
}



void getPredicted() {
  
  TString f_NoMETCut("myHist_1171456__OS_noDupWt_isoDil08_preDil08noIso_hltMu9E15.root");
  TString f_wMETCut("myHist_1433600__OS_noDupWt_isoDil08_preDil08noIso_preMet08_hltMu9E15.root");

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

    
    //Double_t nDYErr2 = p_nllIn_ee.first + pow(p_k.first,2)*p_nIn_emu.first + pow(p_nIn_emu.first,2)*kErr2_ee;
    //Double_t dlEvans_err = pow(predicted_ee,2)*(rOutInErr2_ee/pow(p_rOutIn_ee.first,2) + nDYErr2/pow(p_nllIn_ee.first,2));
    
    //cout << sqrt(dlEvans_err) << endl;
    
  }
  
  cout << endl << "| *nJets*  |  *k*  |" << endl;
  for(unsigned int i = 0 ; i < v_k_ee.size(); i++) 
    cout << "|  " << i << "  |  " 
	 << v_k_ee.at(i).first << " +/- " << sqrt(v_k_ee.at(i).second) << "  |" << endl;
  
  
  cout << endl << endl << "MM" << endl;
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
    

    
    
  

  

#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TIterator.h"
#include <iostream>
#include "TStyle.h"
#include "TChain.h"
#include "TCut.h"

using namespace std;

void makeWWFakeRates(bool doels = true, bool domus = true);

// Method by pointer
TH2F* eff2(TH2F* h1, TH2F* h2, const char* name="eff"){

  // first, verify that all histograms have same binning
  // nx is the number of visible bins
  // nxtot = nx+2 includes underflow and overflow
  Int_t nx = h1->GetNbinsX();
  if (h2->GetNbinsX() != nx) {
    cout << "Histograms must have same number of bins" << endl;
    return 0;
  }

  // get the new histogram
  TH2F* temp = (TH2F*) h1->Clone(name);
  temp->SetTitle(name);
  temp->Reset();
  temp->Sumw2();

  // Do the calculation
  temp->Divide(h2,h1,1.,1.,"B");

  // Done
  return temp;
}


// Method by name
TH2F* eff2(const char* name1, const char* name2, const char* name="eff"){

  // Get a list of object and their iterator
  TList* list = gDirectory->GetList() ;
  TIterator* iter = list->MakeIterator();

  // Loop over objects, set the pointers
  TObject* obj;
  TH2F* h1=0;
  TH2F* h2=0;
  TString str1 = Form("%s",name1);
  TString str2 = Form("%s",name2);
  while((obj=iter->Next())) {
    TString objName = obj->GetName();
    if (objName == str1) h1 = (TH2F*) obj;
    if (objName == str2) h2 = (TH2F*) obj;
  }

  // quit if not found
  if (h1 == 0) {
    cout << "Histogram " << name1 << " not found" << endl;
    return 0;
  }
  if (h2 == 0) {
    cout << "Histogram " << name2 << " not found" << endl;
    return 0;
  }

  // Call the method by pointer
  TH2F* temp = eff2(h1, h2, name);
  return temp;
}

void makeWWFakeRates(bool doels, bool domus) {
    gStyle->SetOptStat(0);

    //--------------------------------------------------------------------------
    // Load lepton data that you want to use to
    // make your fake rate. It should be a baby
    // ntuple.  
    //--------------------------------------------------------------------------

    TChain *ch_el = new TChain("tree");
    TChain *ch_egmon = new TChain("tree");
    TFile* fout_el(0);
    if (doels) {
        ch_el   ->Add("/tas/cms2/FRBabies/FakeRates10November2010-v2/EG.root");
        ch_egmon->Add("/tas/cms2/FRBabies/FakeRates10November2010-v2/EGMon.root");
        fout_el = TFile::Open("ww_el_fr_EGandEGMon.root","RECREATE");
        //ch_el->Add("/tas/cms2/FRBabies/FakeRates10November2010-v2/qcd_pt_30to50_fall10.root");
        //fout_el = TFile::Open("ww_el_fr_qcd_pt_30to50_fall10.root","RECREATE");
    }

    TChain *ch_mu = new TChain("tree");
    TFile* fout_mu(0);
    if (domus) {
        ch_mu->Add("/tas/cms2/FRBabies/FakeRates10November2010-v2/Mu.root");
        fout_mu = TFile::Open("ww_mu_fr_Mu.root","RECREATE");
        //ch_mu->Add("/tas/cms2/FRBabies/FakeRates10November2010-v2/mu15.root");
        //ch_mu->Add("/tas/cms2/FRBabies/FakeRates10November2010-v2/mu15_1.root");
        //fout_mu = TFile::Open("ww_mu_fr_mu15.root","RECREATE");
        //ch_mu->Add("/tas/cms2/FRBabies/FakeRates10November2010-v2/mu10.root");
        //fout_mu = TFile::Open("ww_mu_fr_mu10.root","RECREATE");
    }

    //--------------------------------------------------------------------------
    // Define numerator and denominator cuts
    //--------------------------------------------------------------------------

    // A cut against Ws
    TCut notWCut  = "tcmet<20 && mt<25"; 

    // A pt cut...
    // Remember, we use 35 (at least for muons) to
    // minimize the impact of Ws
    TCut ptCut    = "pt>10 && pt<35";

    // Only consider events with >0 pfjet with pt
    // above a threshold separated by at least dR
    TCut jetCut   = "ptpfj1>15";

    //
    // The trigger selections
    //

    // For EG
    TCut trgCutEl_eg = "((el10_lw>1 || el10_sw>1 || el15_lw>1) && run < 141956) ||\
                        ((el15_sw>1 || el20_sw>1) && run <= 145761)";
    // For EGMon
    TCut trgCutEl_egmon = "((el10_sw>1 || el17_sw>1) && run < 149181) ||\
                           (el10_sw_v2>1 || el17_sw_v2>1)";
    // For Mu
    TCut trgCutMu = "mu9>1 || mu11>1 || mu15>1";

    // Numerator selections
    TCut is_el_num_wwV1 = "num_wwV1&&abs(id)==11"+jetCut+ptCut+notWCut;
    TCut is_mu_num_wwV1 = "num_wwV1&&abs(id)==13"+trgCutMu+jetCut+ptCut+notWCut;

    // Denominator selections
    TCut is_el_v1_wwV1    = "v1_wwV1&&abs(id)==11"+jetCut+ptCut+notWCut;
    TCut is_el_v2_wwV1    = "v2_wwV1&&abs(id)==11"+jetCut+ptCut+notWCut;
    TCut is_el_v3_wwV1    = "v3_wwV1&&abs(id)==11"+jetCut+ptCut+notWCut;
    TCut is_el_v4_wwV1    = "v4_wwV1&&abs(id)==11"+jetCut+ptCut+notWCut;
    TCut is_mu_fo_wwV1_04 = "fo_wwV1_04&&abs(id)==13"+trgCutMu+jetCut+ptCut+notWCut;
    TCut is_mu_fo_wwV1_10 = "fo_wwV1_10&&abs(id)==13"+trgCutMu+jetCut+ptCut+notWCut;
    TCut is_mu_fo_wwV1_10_d0 = "fo_wwV1_10_d0&&abs(id)==13"+trgCutMu+jetCut+ptCut+notWCut;

    //--------------------------------------------------------------------------
    // Define pt and eta bins of fake rate histograms
    //--------------------------------------------------------------------------

    double ybin[6]={10.,15.,20.,25.,30.,35.};
    int nbinsy = 5;
    //double ybin[4]={20.,25.,30.,35.};
    //int nbinsy = 3;
    double xbin[5]={0.0,1.0,1.479,2.0,2.5};
    int nbinsx = 4;

    //--------------------------------------------------------------------------
    // Book numerator and denominator histograms
    //--------------------------------------------------------------------------

    if (doels) fout_el->cd();
    TH2F* el_num_wwV1   = new TH2F("el_num_wwV1","el_num_wwV1",nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v1_wwV1    = new TH2F("el_v1_wwV1","el_v1_wwV1",nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v2_wwV1    = new TH2F("el_v2_wwV1","el_v2_wwV1",nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v3_wwV1    = new TH2F("el_v3_wwV1","el_v3_wwV1",nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v4_wwV1    = new TH2F("el_v4_wwV1","el_v4_wwV1",nbinsx,xbin,nbinsy,ybin);

    if (domus) fout_mu->cd();
    TH2F* mu_num_wwV1   = new TH2F("mu_num_wwV1","mu_num_wwV1",nbinsx,xbin,nbinsy,ybin);
    TH2F* mu_fo_wwV1_04 = new TH2F("mu_fo_wwV1_04","mu_fo_wwV1_04",nbinsx,xbin,nbinsy,ybin);
    TH2F* mu_fo_wwV1_10 = new TH2F("mu_fo_wwV1_10","mu_fo_wwV1_10",nbinsx,xbin,nbinsy,ybin);
    TH2F* mu_fo_wwV1_10_d0 = new TH2F("mu_fo_wwV1_10_d0","mu_fo_wwV1_10_d0",nbinsx,xbin,nbinsy,ybin);

    //--------------------------------------------------------------------------
    // Fill Histograms
    //--------------------------------------------------------------------------

    if (doels) {
        fout_el->cd();
        ch_el->Draw("pt:abs(eta)>>el_num_wwV1",trgCutEl_eg&&is_el_num_wwV1);
        ch_el->Draw("pt:abs(eta)>>el_v1_wwV1",trgCutEl_eg&&is_el_v1_wwV1);
        ch_el->Draw("pt:abs(eta)>>el_v2_wwV1",trgCutEl_eg&&is_el_v2_wwV1);
        ch_el->Draw("pt:abs(eta)>>el_v3_wwV1",trgCutEl_eg&&is_el_v3_wwV1);
        ch_el->Draw("pt:abs(eta)>>el_v4_wwV1",trgCutEl_eg&&is_el_v4_wwV1);
        ch_egmon->Draw("pt:abs(eta)>>+el_num_wwV1",trgCutEl_egmon&&is_el_num_wwV1);
        ch_egmon->Draw("pt:abs(eta)>>+el_v1_wwV1",trgCutEl_egmon&&is_el_v1_wwV1);
        ch_egmon->Draw("pt:abs(eta)>>+el_v2_wwV1",trgCutEl_egmon&&is_el_v2_wwV1);
        ch_egmon->Draw("pt:abs(eta)>>+el_v3_wwV1",trgCutEl_egmon&&is_el_v3_wwV1);
        ch_egmon->Draw("pt:abs(eta)>>+el_v4_wwV1",trgCutEl_egmon&&is_el_v4_wwV1);
    }

    if (domus) {
        fout_mu->cd();
        ch_mu->Draw("pt:abs(eta)>>mu_num_wwV1",is_mu_num_wwV1);
        ch_mu->Draw("pt:abs(eta)>>mu_fo_wwV1_04",is_mu_fo_wwV1_04);
        ch_mu->Draw("pt:abs(eta)>>mu_fo_wwV1_10",is_mu_fo_wwV1_10);
        ch_mu->Draw("pt:abs(eta)>>mu_fo_wwV1_10_d0",is_mu_fo_wwV1_10_d0);
    }

    //--------------------------------------------------------------------------
    // Make fake rate histograms. Last argument is
    // the histogram name
    //--------------------------------------------------------------------------

    if (doels){
      fout_el->cd();
      eff2(el_v1_wwV1,el_num_wwV1,"el_fr_v1_wwV1");
      eff2(el_v2_wwV1,el_num_wwV1,"el_fr_v2_wwV1");
      eff2(el_v3_wwV1,el_num_wwV1,"el_fr_v3_wwV1");
      eff2(el_v4_wwV1,el_num_wwV1,"el_fr_v4_wwV1");
      fout_el->Write();
      fout_el->Close();
    }

    if (domus){
      fout_mu->cd();
      eff2(mu_fo_wwV1_04,mu_num_wwV1,"mu_fr_fo_wwV1_04");
      eff2(mu_fo_wwV1_10,mu_num_wwV1,"mu_fr_fo_wwV1_10");
      eff2(mu_fo_wwV1_10_d0,mu_num_wwV1,"mu_fr_fo_wwV1_10_d0");
      fout_mu->Write();
      fout_mu->Close();
    }
}

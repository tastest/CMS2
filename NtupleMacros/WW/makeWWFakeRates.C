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
    TChain *ch_el = new TChain("tree");
    TFile* fout_el(0);
    if (doels) {
      ch_el   ->Add("/smurf/dmytro/samples/fake_ntuples/FakeRates18May2011/DoubleElectron.root");
      fout_el = TFile::Open("ww_el_fr_smurfV4.root","RECREATE");
    }

    TChain *ch_mu = new TChain("tree");
    TFile* fout_mu(0);
    if (domus) {
        ch_mu->Add("/smurf/dmytro/samples/fake_ntuples/FakeRates18May2011/SingleMu.root");
        fout_mu = TFile::Open("ww_mu_fr_smurfV4.root","RECREATE");
    }
    // A cut against Ws
    TCut notWCut  = "tcmet<20 && mt<25"; 

    // A pt cut...
    // Remember, we use 35 (at least for muons) to
    // minimize the impact of Ws
    TCut ptCut    = "pt>10 && pt<35";

    // Only consider events with >0 pfjet with pt
    // above a threshold separated by at least dR
    TCut jetCut5    = "ptpfj1>5";
    TCut jetCut15   = "ptpfj1>15";
    TCut jetCut30   = "ptpfj1>30";
    TCut jetCut50   = "ptpfj1>50";

    //
    // The trigger selections
    //
    TCut trgCutEl = "ele8_vstar>1 || ele8_CaloIdL_CaloIsoVL_vstar>1 || ele17_CaloIdL_CaloIsoVL_vstar>1";
    TCut trgCutMu = "mu5_vstar>1 || mu8_vstar>1 || mu12_vstar>1 || mu15_vstar>1";

    // Numerator selections
    TCut is_el_num = "num_el_smurfV3&&abs(id)==11"+trgCutEl+ptCut+notWCut;
    TCut is_mu_num = "num_mu_smurfV3&&abs(id)==13"+trgCutMu+ptCut+notWCut;

    // Denominator selections
    TCut is_el_fo_v1 = "v1_el_smurfV1 && abs(id)==11" +trgCutEl+ptCut+notWCut;
    TCut is_el_fo_v2 = "v2_el_smurfV1 && abs(id)==11" +trgCutEl+ptCut+notWCut;
    TCut is_el_fo_v3 = "v3_el_smurfV1 && abs(id)==11" +trgCutEl+ptCut+notWCut;
    TCut is_el_fo_v4 = "v4_el_smurfV1 && abs(id)==11" +trgCutEl+ptCut+notWCut;
    TCut is_mu_fo_m1 = "fo_mu_smurf_10 && abs(id)==13"+trgCutMu+ptCut+notWCut;
    TCut is_mu_fo_m2 = "fo_mu_smurf_04 && abs(id)==13"+trgCutMu+ptCut+notWCut;

    //--------------------------------------------------------------------------
    // Define pt and eta bins of fake rate histograms
    //--------------------------------------------------------------------------

    double ybin[6]={10.,15.,20.,25.,30.,35.};
    int nbinsy = 5;
    double xbin[5]={0.0,1.0,1.479,2.0,2.5};
    int nbinsx = 4;

    //--------------------------------------------------------------------------
    // Book numerator and denominator histograms
    //--------------------------------------------------------------------------

    if (doels) fout_el->cd();
    TH2F* el_num_5    = new TH2F("el_num_5", "Electrons passed final selection (JetPt>5)", nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v1_5     = new TH2F("el_v1_5",  "Electrons passed V1 selection (JetPt>5)",    nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v2_5     = new TH2F("el_v2_5",  "Electrons passed V2 selection (JetPt>5)",    nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v3_5     = new TH2F("el_v3_5",  "Electrons passed V3 selection (JetPt>5)",    nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v4_5     = new TH2F("el_v4_5",  "Electrons passed V4 selection (JetPt>5)",    nbinsx,xbin,nbinsy,ybin);
    TH2F* el_num_15   = new TH2F("el_num_15","Electrons passed final selection (JetPt>15)", nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v1_15    = new TH2F("el_v1_15", "Electrons passed V1 selection (JetPt>15)",    nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v2_15    = new TH2F("el_v2_15", "Electrons passed V2 selection (JetPt>15)",    nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v3_15    = new TH2F("el_v3_15", "Electrons passed V3 selection (JetPt>15)",    nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v4_15    = new TH2F("el_v4_15", "Electrons passed V4 selection (JetPt>15)",    nbinsx,xbin,nbinsy,ybin);
    TH2F* el_num_30   = new TH2F("el_num_30","Electrons passed final selection (JetPt>30)", nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v1_30    = new TH2F("el_v1_30", "Electrons passed V1 selection (JetPt>30)",    nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v2_30    = new TH2F("el_v2_30", "Electrons passed V2 selection (JetPt>30)",    nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v3_30    = new TH2F("el_v3_30", "Electrons passed V3 selection (JetPt>30)",    nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v4_30    = new TH2F("el_v4_30", "Electrons passed V4 selection (JetPt>30)",    nbinsx,xbin,nbinsy,ybin);
    TH2F* el_num_50   = new TH2F("el_num_50","Electrons passed final selection (JetPt>50)", nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v1_50    = new TH2F("el_v1_50", "Electrons passed V1 selection (JetPt>50)",    nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v2_50    = new TH2F("el_v2_50", "Electrons passed V2 selection (JetPt>50)",    nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v3_50    = new TH2F("el_v3_50", "Electrons passed V3 selection (JetPt>50)",    nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v4_50    = new TH2F("el_v4_50", "Electrons passed V4 selection (JetPt>50)",    nbinsx,xbin,nbinsy,ybin);

    if (domus) fout_mu->cd();
    TH2F* mu_num_5    = new TH2F("mu_num_5", "Muons passed final selection (JetPt>5)",     nbinsx,xbin,nbinsy,ybin);
    TH2F* mu_m1_5     = new TH2F("mu_m1_5",  "Muons passed M1 selection (JetPt>5)",        nbinsx,xbin,nbinsy,ybin);
    TH2F* mu_m2_5     = new TH2F("mu_m2_5",  "Muons passed M2 selection (JetPt>5)",        nbinsx,xbin,nbinsy,ybin);
    TH2F* mu_num_15   = new TH2F("mu_num_15","Muons passed final selection (JetPt>15)",     nbinsx,xbin,nbinsy,ybin);
    TH2F* mu_m1_15    = new TH2F("mu_m1_15", "Muons passed M1 selection (JetPt>15)",        nbinsx,xbin,nbinsy,ybin);
    TH2F* mu_m2_15    = new TH2F("mu_m2_15", "Muons passed M2 selection (JetPt>15)",        nbinsx,xbin,nbinsy,ybin);
    TH2F* mu_num_30   = new TH2F("mu_num_30","Muons passed final selection (JetPt>30)",     nbinsx,xbin,nbinsy,ybin);
    TH2F* mu_m1_30    = new TH2F("mu_m1_30", "Muons passed M1 selection (JetPt>30)",        nbinsx,xbin,nbinsy,ybin);
    TH2F* mu_m2_30    = new TH2F("mu_m2_30", "Muons passed M2 selection (JetPt>30)",        nbinsx,xbin,nbinsy,ybin);
    TH2F* mu_num_50   = new TH2F("mu_num_50","Muons passed final selection (JetPt>50)",     nbinsx,xbin,nbinsy,ybin);
    TH2F* mu_m1_50    = new TH2F("mu_m1_50", "Muons passed M1 selection (JetPt>50)",        nbinsx,xbin,nbinsy,ybin);
    TH2F* mu_m2_50    = new TH2F("mu_m2_50", "Muons passed M2 selection (JetPt>50)",        nbinsx,xbin,nbinsy,ybin);

    //--------------------------------------------------------------------------
    // Fill Histograms
    //--------------------------------------------------------------------------

    if (doels) {
        fout_el->cd();
        ch_el->Draw("pt:abs(eta)>>el_num_5", is_el_num+jetCut5);
        ch_el->Draw("pt:abs(eta)>>el_v1_5",  is_el_fo_v1+jetCut5);
        ch_el->Draw("pt:abs(eta)>>el_v2_5",  is_el_fo_v2+jetCut5);
        ch_el->Draw("pt:abs(eta)>>el_v3_5",  is_el_fo_v3+jetCut5);
        ch_el->Draw("pt:abs(eta)>>el_v4_5",  is_el_fo_v4+jetCut5);
        ch_el->Draw("pt:abs(eta)>>el_num_15", is_el_num+jetCut15);
        ch_el->Draw("pt:abs(eta)>>el_v1_15",  is_el_fo_v1+jetCut15);
        ch_el->Draw("pt:abs(eta)>>el_v2_15",  is_el_fo_v2+jetCut15);
        ch_el->Draw("pt:abs(eta)>>el_v3_15",  is_el_fo_v3+jetCut15);
        ch_el->Draw("pt:abs(eta)>>el_v4_15",  is_el_fo_v4+jetCut15);
        ch_el->Draw("pt:abs(eta)>>el_num_30", is_el_num+jetCut30);
        ch_el->Draw("pt:abs(eta)>>el_v1_30",  is_el_fo_v1+jetCut30);
        ch_el->Draw("pt:abs(eta)>>el_v2_30",  is_el_fo_v2+jetCut30);
        ch_el->Draw("pt:abs(eta)>>el_v3_30",  is_el_fo_v3+jetCut30);
        ch_el->Draw("pt:abs(eta)>>el_v4_30",  is_el_fo_v4+jetCut30);
        ch_el->Draw("pt:abs(eta)>>el_num_50", is_el_num+jetCut50);
        ch_el->Draw("pt:abs(eta)>>el_v1_50",  is_el_fo_v1+jetCut50);
        ch_el->Draw("pt:abs(eta)>>el_v2_50",  is_el_fo_v2+jetCut50);
        ch_el->Draw("pt:abs(eta)>>el_v3_50",  is_el_fo_v3+jetCut50);
        ch_el->Draw("pt:abs(eta)>>el_v4_50",  is_el_fo_v4+jetCut50);
    }

    if (domus) {
        fout_mu->cd();
        ch_mu->Draw("pt:abs(eta)>>mu_num_5", is_mu_num+jetCut5);
        ch_mu->Draw("pt:abs(eta)>>mu_m1_5",  is_mu_fo_m1+jetCut5);
        ch_mu->Draw("pt:abs(eta)>>mu_m2_5",  is_mu_fo_m2+jetCut5);
        ch_mu->Draw("pt:abs(eta)>>mu_num_15", is_mu_num+jetCut15);
        ch_mu->Draw("pt:abs(eta)>>mu_m1_15",  is_mu_fo_m1+jetCut15);
        ch_mu->Draw("pt:abs(eta)>>mu_m2_15",  is_mu_fo_m2+jetCut15);
        ch_mu->Draw("pt:abs(eta)>>mu_num_30", is_mu_num+jetCut30);
        ch_mu->Draw("pt:abs(eta)>>mu_m1_30",  is_mu_fo_m1+jetCut30);
        ch_mu->Draw("pt:abs(eta)>>mu_m2_30",  is_mu_fo_m2+jetCut30);
        ch_mu->Draw("pt:abs(eta)>>mu_num_50", is_mu_num+jetCut50);
        ch_mu->Draw("pt:abs(eta)>>mu_m1_50",  is_mu_fo_m1+jetCut50);
        ch_mu->Draw("pt:abs(eta)>>mu_m2_50",  is_mu_fo_m2+jetCut50);
    }

    //--------------------------------------------------------------------------
    // Make fake rate histograms. Last argument is
    // the histogram name
    //--------------------------------------------------------------------------

    if (doels){
      fout_el->cd();
      eff2(el_v1_5,  el_num_5,  "el_fr_v1_5");
      eff2(el_v2_5,  el_num_5,  "el_fr_v2_5");
      eff2(el_v3_5,  el_num_5,  "el_fr_v3_5");
      eff2(el_v4_5,  el_num_5,  "el_fr_v4_5");
      eff2(el_v1_15, el_num_15, "el_fr_v1_15");
      eff2(el_v2_15, el_num_15, "el_fr_v2_15");
      eff2(el_v3_15, el_num_15, "el_fr_v3_15");
      eff2(el_v4_15, el_num_15, "el_fr_v4_15");
      eff2(el_v1_30, el_num_30, "el_fr_v1_30");
      eff2(el_v2_30, el_num_30, "el_fr_v2_30");
      eff2(el_v3_30, el_num_30, "el_fr_v3_30");
      eff2(el_v4_30, el_num_30, "el_fr_v4_30");
      eff2(el_v1_50, el_num_50, "el_fr_v1_50");
      eff2(el_v2_50, el_num_50, "el_fr_v2_50");
      eff2(el_v3_50, el_num_50, "el_fr_v3_50");
      eff2(el_v4_50, el_num_50, "el_fr_v4_50");
      fout_el->Write();
      fout_el->Close();
    }

    if (domus){
      fout_mu->cd();
      eff2(mu_m1_5,  mu_num_5,  "mu_fr_m1_5");
      eff2(mu_m2_5,  mu_num_5,  "mu_fr_m2_5");
      eff2(mu_m1_15, mu_num_15, "mu_fr_m1_15");
      eff2(mu_m2_15, mu_num_15, "mu_fr_m2_15");
      eff2(mu_m1_30, mu_num_30, "mu_fr_m1_30");
      eff2(mu_m2_30, mu_num_30, "mu_fr_m2_30");
      eff2(mu_m1_50, mu_num_50, "mu_fr_m1_50");
      eff2(mu_m2_50, mu_num_50, "mu_fr_m2_50");
      fout_mu->Write();
      fout_mu->Close();
    }
}

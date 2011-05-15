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
      ch_el   ->Add("/smurf/dmytro/samples/fake_ntuples/FakeRates8May2011/DoubleElectron.root");
      fout_el = TFile::Open("ww_el_fr.root","RECREATE");
    }

    TChain *ch_mu = new TChain("tree");
    TFile* fout_mu(0);
    if (domus) {
        ch_mu->Add("/smurf/dmytro/samples/fake_ntuples/FakeRates8May2011/SingleMu.root");
        fout_mu = TFile::Open("ww_mu_fr.root","RECREATE");
    }
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
    TCut trgCutEl = "ele8_vstar>1 || ele8_CaloIdL_CaloIsoVL_vstar>1 || ele17_CaloIdL_CaloIsoVL_vstar>1";
    TCut trgCutMu = "mu5_vstar>1 || mu8_vstar>1 || mu12_vstar>1 || mu15_vstar>1";

    // Numerator selections
    TCut is_el_num = "num_el_smurfV3&&abs(id)==11"+trgCutEl+jetCut+ptCut+notWCut;
    TCut is_mu_num = "num_mu_smurfV3&&abs(id)==13"+trgCutMu+jetCut+ptCut+notWCut;

    // Denominator selections
    TCut is_el_fo_v1 = "v1_el_smurfV1 && abs(id)==11" +trgCutEl+jetCut+ptCut+notWCut;
    TCut is_el_fo_v2 = "v2_el_smurfV1 && abs(id)==11" +trgCutEl+jetCut+ptCut+notWCut;
    TCut is_el_fo_v3 = "v3_el_smurfV1 && abs(id)==11" +trgCutEl+jetCut+ptCut+notWCut;
    TCut is_el_fo_v4 = "v4_el_smurfV1 && abs(id)==11" +trgCutEl+jetCut+ptCut+notWCut;
    TCut is_mu_fo_m1 = "fo_mu_smurf_10 && abs(id)==13"+trgCutMu+jetCut+ptCut+notWCut;
    TCut is_mu_fo_m2 = "fo_mu_smurf_04 && abs(id)==13"+trgCutMu+jetCut+ptCut+notWCut;

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
    TH2F* el_num   = new TH2F("el_num","Electrons passed final selection", nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v1    = new TH2F("el_v1", "Electrons passed V1 selection",    nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v2    = new TH2F("el_v2", "Electrons passed V2 selection",    nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v3    = new TH2F("el_v3", "Electrons passed V3 selection",    nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v4    = new TH2F("el_v4", "Electrons passed V4 selection",    nbinsx,xbin,nbinsy,ybin);

    if (domus) fout_mu->cd();
    TH2F* mu_num   = new TH2F("mu_num","Muons passed final selection",     nbinsx,xbin,nbinsy,ybin);
    TH2F* mu_m1    = new TH2F("mu_m1", "Muons passed M1 selection",        nbinsx,xbin,nbinsy,ybin);
    TH2F* mu_m2    = new TH2F("mu_m2", "Muons passed M2 selection",        nbinsx,xbin,nbinsy,ybin);

    //--------------------------------------------------------------------------
    // Fill Histograms
    //--------------------------------------------------------------------------

    if (doels) {
        fout_el->cd();
        ch_el->Draw("pt:abs(eta)>>el_num", is_el_num);
        ch_el->Draw("pt:abs(eta)>>el_v1",  is_el_fo_v1);
        ch_el->Draw("pt:abs(eta)>>el_v2",  is_el_fo_v2);
        ch_el->Draw("pt:abs(eta)>>el_v3",  is_el_fo_v3);
        ch_el->Draw("pt:abs(eta)>>el_v4",  is_el_fo_v4);
    }

    if (domus) {
        fout_mu->cd();
        ch_mu->Draw("pt:abs(eta)>>mu_num", is_mu_num);
        ch_mu->Draw("pt:abs(eta)>>mu_m1",  is_mu_fo_m1);
        ch_mu->Draw("pt:abs(eta)>>mu_m2",  is_mu_fo_m2);
    }

    //--------------------------------------------------------------------------
    // Make fake rate histograms. Last argument is
    // the histogram name
    //--------------------------------------------------------------------------

    if (doels){
      fout_el->cd();
      eff2(el_v1, el_num, "el_fr_v1");
      eff2(el_v2, el_num, "el_fr_v2");
      eff2(el_v3, el_num, "el_fr_v3");
      eff2(el_v4, el_num, "el_fr_v4");
      fout_el->Write();
      fout_el->Close();
    }

    if (domus){
      fout_mu->cd();
      eff2(mu_m1, mu_num, "mu_fr_m1");
      eff2(mu_m2, mu_num, "mu_fr_m2");
      fout_mu->Write();
      fout_mu->Close();
    }
}

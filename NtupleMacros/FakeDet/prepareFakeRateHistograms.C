#include "TFile.h"
#include "TH2F.h"

const char* makeName(const char* bin, const char* name) {
  return Form("%s_h%s",bin,name);
}

void prepareFakeRateHistograms() {

  vector<char*> qcdBins;
  vector<float> qcdBinXSec;
  vector<float> qcdBinFilterEff;
  vector<float> qcdBinNEvents; // number of processed events corresponding to used ntuples

  vector<char*> histograms;

  qcdBins.push_back("qcd_0_15");
  qcdBins.push_back("qcd_15_20");
  qcdBins.push_back("qcd_20_30");
  qcdBins.push_back("qcd_30_50");
  qcdBins.push_back("qcd_50_80");
  qcdBins.push_back("qcd_80_120");
  qcdBins.push_back("qcd_120_170");
  qcdBins.push_back("qcd_170_230");
  qcdBins.push_back("qcd_230_300");
  qcdBins.push_back("qcd_300_380");
  qcdBins.push_back("qcd_380_470");
  qcdBins.push_back("qcd_470_600");

  qcdBinFilterEff.push_back(0.964);
  qcdBinFilterEff.push_back(1.);
  qcdBinFilterEff.push_back(1.);
  qcdBinFilterEff.push_back(1.);
  qcdBinFilterEff.push_back(1.);
  qcdBinFilterEff.push_back(1.);
  qcdBinFilterEff.push_back(1.);
  qcdBinFilterEff.push_back(1.);
  qcdBinFilterEff.push_back(1.);
  qcdBinFilterEff.push_back(1.);
  qcdBinFilterEff.push_back(1.);
  qcdBinFilterEff.push_back(1.);

  qcdBinXSec.push_back(55000000000.00);
  qcdBinXSec.push_back(1460000000.00);
  qcdBinXSec.push_back(630000000.00);
  qcdBinXSec.push_back(163000000.00);
  qcdBinXSec.push_back(21600000.00);
  qcdBinXSec.push_back(3080000.00);
  qcdBinXSec.push_back(494000.00);
  qcdBinXSec.push_back(101000.00);
  qcdBinXSec.push_back(24500.00);
  qcdBinXSec.push_back(6240.00);
  qcdBinXSec.push_back(1780.00);
  qcdBinXSec.push_back(683.);


  qcdBinNEvents.push_back(611787);
  qcdBinNEvents.push_back(1286976);
  qcdBinNEvents.push_back(1908861);
  qcdBinNEvents.push_back(1160479);
  qcdBinNEvents.push_back(914740);
  qcdBinNEvents.push_back(1258762);
  qcdBinNEvents.push_back(1260951);
  qcdBinNEvents.push_back(934870);
  qcdBinNEvents.push_back(799844);
  qcdBinNEvents.push_back(1274039);
  qcdBinNEvents.push_back(1246217);
  qcdBinNEvents.push_back(1315614);

  histograms.push_back("den_ele");
  histograms.push_back("den_wo_leading_ele");
  histograms.push_back("den_wo_second_leading_ele");

  histograms.push_back("num_ell");
  histograms.push_back("num_wo_leading_ell");
  histograms.push_back("num_wo_second_leading_ell");

  histograms.push_back("num_elt");
  histograms.push_back("num_wo_leading_elt");
  histograms.push_back("num_wo_second_leading_elt");

  TFile *file = TFile::Open("fakeRatesFull.root","UPDATE");

  for ( unsigned int histo = 0;
	histo < histograms.size();
	++histo ) {
    bool first = true;
    TH2F* clone = 0;
    for ( unsigned int bin = 0;
	  bin < qcdBins.size();
	  ++bin ) {

      TH2F *tmp = dynamic_cast<TH2F*>(file->Get(makeName(qcdBins[bin],histograms[histo])));
      if ( tmp == 0. ) continue;
      tmp->SetDirectory(0);
      if ( first ) {
	clone = dynamic_cast<TH2F*>tmp->Clone(histograms[histo]);
	clone->Scale(qcdBinXSec[bin]*qcdBinFilterEff[bin]/qcdBinNEvents[bin]);
	first = false;
      } else {
	clone->Add(tmp,qcdBinXSec[bin]*qcdBinFilterEff[bin]/qcdBinNEvents[bin]);
      }
    }
  }

  // create fake rate histograms, custom-made
  TH2F* den0  = dynamic_cast<TH2F*>(file->Get(histograms[0]));
  TH2F* den1  = dynamic_cast<TH2F*>(file->Get(histograms[1]));
  TH2F* den2  = dynamic_cast<TH2F*>(file->Get(histograms[2]));

  TH2F* num0 = dynamic_cast<TH2F*>(file->Get(histograms[3]));
  TH2F* fake0 = dynamic_cast<TH2F*>(num0->Clone(TString(histograms[3]).ReplaceAll("num_","fakeRate_")));
  fake0->Reset();
  fake0->Divide(num0,den0,1.,1.,"B");

  TH2F* num1 = dynamic_cast<TH2F*>(file->Get(histograms[4]));
  TH2F* fake1 = dynamic_cast<TH2F*>(num1->Clone(TString(histograms[4]).ReplaceAll("num_","fakeRate_")));
  fake1->Reset();
  fake1->Divide(num1,den1,1.,1.,"B");

  TH2F* num2 = dynamic_cast<TH2F*>(file->Get(histograms[5]));
  TH2F* fake2 = dynamic_cast<TH2F*>(num2->Clone(TString(histograms[5]).ReplaceAll("num_","fakeRate_")));
  fake2->Reset();
  fake2->Divide(num2,den2,1.,1.,"B");

  TH2F* num3 = dynamic_cast<TH2F*>(file->Get(histograms[6]));
  TH2F* fake3 = dynamic_cast<TH2F*>(num3->Clone(TString(histograms[6]).ReplaceAll("num_","fakeRate_")));
  fake3->Reset();
  fake3->Divide(num3,den0,1.,1.,"B");

  TH2F* num4 = dynamic_cast<TH2F*>(file->Get(histograms[7]));
  TH2F* fake4 = dynamic_cast<TH2F*>(num4->Clone(TString(histograms[7]).ReplaceAll("num_","fakeRate_")));
  fake4->Reset();
  fake4->Divide(num4,den1,1.,1.,"B");

  TH2F* num5 = dynamic_cast<TH2F*>(file->Get(histograms[8]));
  TH2F* fake5 = dynamic_cast<TH2F*>(num5->Clone(TString(histograms[8]).ReplaceAll("num_","fakeRate_")));
  fake5->Reset();
  fake5->Divide(num5,den2,1.,1.,"B");

  file->Write();
  file->Close();

}

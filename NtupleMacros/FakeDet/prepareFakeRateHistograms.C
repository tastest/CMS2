#include "TFile.h"
#include "TH2F.h"
#include <vector>

using std::vector;

const char* makeName(const char* bin, const char* name) {
  return Form("%s_h%s",bin,name);
}

void prepareFakeRateHistograms() {

  vector<char*> qcdBins;
  vector<char*> histograms;

  qcdBins.push_back("QCDEMenrichedPt20to30");
  qcdBins.push_back("QCDEMenrichedPt30to80");
  qcdBins.push_back("QCDEMenrichedPt80to170");

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
	first = false;
      } else {
	clone->Add(tmp,1.);
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

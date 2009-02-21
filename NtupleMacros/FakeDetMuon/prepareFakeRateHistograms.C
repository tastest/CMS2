#include "TFile.h"
#include "TH2F.h"
#include <vector>
#include "checkZeroBin2D.C"

using std::vector;

const char* makeName(const char* bin, const char* name) {
  return Form("%s_h%s",bin,name);
}

void prepareFakeRateHistograms() {

  vector<char*> qcdBins;
  vector<char*> histograms;

  bool useQCDBCtoE = true;
  bool useInclusiveSamples = false;
  bool useInclusiveMuonSamples = true;

  if(useInclusiveSamples && useInclusiveMuonSamples) {
    cout<<"Both incl QCD and inclusiveMuon samples are set to true - currently not possible. EXIT"<<endl;
    break;
  }
    

  if ( useInclusiveSamples ) {
    qcdBins.push_back("QCDpt30");
  }
  else if( useInclusiveMuonSamples ) {
    qcdBins.push_back("InclusiveMuPt15");
  }
  else {
    qcdBins.push_back("QCDEMenrichedPt20to30");
    qcdBins.push_back("QCDEMenrichedPt30to80");
    qcdBins.push_back("QCDEMenrichedPt80to170");

    if(useQCDBCtoE) {
      qcdBins.push_back("QCDBCtoEPt20to30");
      qcdBins.push_back("QCDBCtoEPt30to80");
      qcdBins.push_back("QCDBCtoEPt80to170");
    }
  }

  histograms.push_back("den_muo");
  histograms.push_back("den_wo_leading_muo");
  histograms.push_back("den_wo_second_leading_muo");

  histograms.push_back("num_mll");
  histograms.push_back("num_wo_leading_mll");
  histograms.push_back("num_wo_second_leading_mll");

  histograms.push_back("num_mlt");
  histograms.push_back("num_wo_leading_mlt");
  histograms.push_back("num_wo_second_leading_mlt");

  TFile *file = TFile::Open("MuoFakes.root","UPDATE");

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
  
  checkZeroBin2D((TH2D*)den0);
  checkZeroBin2D((TH2D*)den1);
  checkZeroBin2D((TH2D*)den2);
  
  TH2F* num0 = dynamic_cast<TH2F*>(file->Get(histograms[3]));
  TH2F* fake0 = dynamic_cast<TH2F*>(num0->Clone(TString(histograms[3]).ReplaceAll("num_","fakeRate_")));

  checkZeroBin2D((TH2D*)num0);
  checkZeroBin2D((TH2D*)fake0);

  fake0->Reset();
  fake0->Divide(num0,den0,1.,1.,"B");

  TH2F* num1 = dynamic_cast<TH2F*>(file->Get(histograms[4]));
  TH2F* fake1 = dynamic_cast<TH2F*>(num1->Clone(TString(histograms[4]).ReplaceAll("num_","fakeRate_")));
  checkZeroBin2D((TH2D*)num1);
  checkZeroBin2D((TH2D*)fake1);
  fake1->Reset();
  fake1->Divide(num1,den1,1.,1.,"B");

  TH2F* num2 = dynamic_cast<TH2F*>(file->Get(histograms[5]));
  TH2F* fake2 = dynamic_cast<TH2F*>(num2->Clone(TString(histograms[5]).ReplaceAll("num_","fakeRate_")));
  checkZeroBin2D((TH2D*)num2);
  checkZeroBin2D((TH2D*)fake2);
  fake2->Reset();
  fake2->Divide(num2,den2,1.,1.,"B");

  TH2F* num3 = dynamic_cast<TH2F*>(file->Get(histograms[6]));
  TH2F* fake3 = dynamic_cast<TH2F*>(num3->Clone(TString(histograms[6]).ReplaceAll("num_","fakeRate_")));
  checkZeroBin2D((TH2D*)num3);
  checkZeroBin2D((TH2D*)fake3);
  fake3->Reset();
  fake3->Divide(num3,den0,1.,1.,"B");

  TH2F* num4 = dynamic_cast<TH2F*>(file->Get(histograms[7]));
  TH2F* fake4 = dynamic_cast<TH2F*>(num4->Clone(TString(histograms[7]).ReplaceAll("num_","fakeRate_")));
  checkZeroBin2D((TH2D*)num4);
  checkZeroBin2D((TH2D*)fake4);
  fake4->Reset();
  fake4->Divide(num4,den1,1.,1.,"B");

  TH2F* num5 = dynamic_cast<TH2F*>(file->Get(histograms[8]));
  TH2F* fake5 = dynamic_cast<TH2F*>(num5->Clone(TString(histograms[8]).ReplaceAll("num_","fakeRate_")));
  fake5->Reset();
  fake5->Divide(num5,den2,1.,1.,"B");

  file->Write();
  file->Close();

}

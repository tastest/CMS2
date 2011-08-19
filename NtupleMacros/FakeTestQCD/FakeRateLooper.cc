#include <sstream>
#include <math.h>
#include "TVector3.h"
#include "TDirectory.h"
#include "CORE/selections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Tools/fakerates.h"
#include "FakeRateLooper.h"

FakeRateLooper::FakeRateLooper (Sample s, cuts_t c, const char *fname) 
  : Looper(s, c, fname)
{
  // zero out the candidate counters (don't comment this out)
  cands_passing_syst_hi_ = 0;
  cands_passing_syst_lo_ = 0;
  cands_passing_event_weight_only_ = 0;
  cands_passing_event_weight_only_w2_ = 0;
}

void FakeRateLooper::BookHistos 	()
{
  Looper::BookHistos();

  gDirectory = histo_directory;

  fake_syst = new TH2F(Form("%s_%s",sample_.name.c_str(),"fake_syst"),Form("%s_%s",sample_.name.c_str(),"fake syst uncertainty;#eta;pt"), 
		       fakeRate().GetNbinsX(), fakeRate().GetXaxis()->GetXbins()->GetArray(),
		       fakeRate().GetNbinsY(), fakeRate().GetYaxis()->GetXbins()->GetArray());

  const unsigned int ptNBins = 16;
  float ptBins[ptNBins+1];
  for ( unsigned int ptBin = 0;
	ptBin <= ptNBins;
	++ptBin) {
    ptBins[ptBin] = float(ptBin)*160./16.;
  }

  const unsigned int etaNBins = 12;
  float etaBins[etaNBins+1];
  for ( unsigned int etaBin = 0;
	etaBin <= etaNBins;
	++etaBin) {
    etaBins[etaBin] = float(etaBin)/2.-3.;
  }

  // determine bin arrays for fakeRate histogram
  const unsigned int fakeXNBins = fakeRate().GetXaxis()->GetNbins();
  float fakeXBins[fakeXNBins+1];
  extractBins(fakeRate().GetXaxis(),fakeXNBins,fakeXBins);
  const unsigned int fakeYNBins = fakeRate().GetYaxis()->GetNbins();
  float fakeYBins[fakeYNBins+1];
  extractBins(fakeRate().GetYaxis(),fakeYNBins,fakeYBins);

  helPt3D_ = book3DVarHist(Form("%s_%s",sample_.name.c_str(),"elPt3D"),
			   Form("%s_%s",sample_.name.c_str(),"elPt3D"),
			   ptNBins,ptBins,
			   fakeXNBins,fakeXBins,
			   fakeYNBins,fakeYBins,
			   "p_{T}^{e} [GeV]","#eta", "p_{T} [GeV]",sample_.histo_color);
  helEta3D_ = book3DVarHist(Form("%s_%s",sample_.name.c_str(),"elEta3D"),
			    Form("%s_%s",sample_.name.c_str(),"elEta3D"),
			    etaNBins,etaBins,
			    fakeXNBins,fakeXBins,
			    fakeYNBins,fakeYBins,
			    "#eta^{e} [GeV]","#eta", "p_{T} [GeV]",sample_.histo_color);
}

cuts_t FakeRateLooper::EventSelect ()
{
  // this doesn't do anything special
  return Looper::EventSelect();
}

void FakeRateLooper::FillEventHistos ()
{
  cuts_t cuts_passed = EventSelect();
     
  if ((cuts_passed & cuts_) == cuts_) {

    // electron loop 
    for ( unsigned int ele = 0; 
          ele < cms2.els_p4().size(); 
          ++ele ) { 
      if ( isFakeable(ele) && !isNumeratorElectron(ele)) {

	const double weight_hi = Weight(ele, 1);
	const double weight_lo = Weight(ele, -1);
	cands_passing_syst_hi_ += weight_hi;
	cands_passing_syst_lo_ += weight_lo;
	const double weight = Looper::Weight(ele);
	cands_passing_event_weight_only_ += weight;
	cands_passing_event_weight_only_w2_ += weight * weight;

	const double err = elFakeProb(ele, 1) - elFakeProb(ele, 0);
	const double eta = cms2.els_p4()[ele].eta();
	const double pt = cms2.els_p4()[ele].pt();
	fake_syst->Fill(eta, pt, weight * err);
	helPt3D_->Fill(cms2.els_p4()[ele].pt(),eta,pt, weight * err);
	helEta3D_->Fill(cms2.els_p4()[ele].eta(),eta,pt, weight * err);
      }
    }
  }

  Looper::FillEventHistos();
}

double FakeRateLooper::Weight (int idx)
{ 
  return Weight(idx, 0);
}

double FakeRateLooper::Weight (int idx, int n_sig_syst)
{
  double weight = Looper::Weight(idx);
  double fr = elFakeProb(idx, n_sig_syst);

  return weight * fr / (1 - fr); 
}

double FakeRateLooper::FakeSyst () const
{
  double err2 = 0;
  for (int i = 0; i <= fake_syst->GetNbinsX() + 1; ++i) {
    for (int j = 0; j <= fake_syst->GetNbinsY() + 1; ++j) {
      double err = fake_syst->GetBinContent(i, j);
      err2 += err * err;
    }
  }
  return sqrt(err2);
}

void FakeRateLooper::End ()
{
  //------------------------------------------------------------
  //Example status message at the end of a looper; edit for your
  //application
  //------------------------------------------------------------

  // treat errors correctly
  fillErrorInPrediction(helPt_,helPt3D_);
  fillErrorInPrediction(helEta_,helEta3D_);

  ostringstream stream;
  
  stream << endl << "=========" << endl;
  stream << "Sample: " << SampleName().c_str() << endl;
  stream << "=========" << endl;
  stream << "Total candidate count: " << CandsCount() << endl;
  stream << "Total weighted events (event weight only):" << fixed << setprecision(1) << CandsPassingEventWeightOnly() << "+-" << RMSEventWeightOnly() << endl;
  stream << "Total weighted events (weight*fr/(1-fr)): ee: " << fixed << setprecision(1) << CandsPassing() << "+-" << RMS() << endl;
  stream << "=========" << endl;
  
  cout << stream.str();

  int ret = fprintf(logfile_, stream.str().c_str());
  if (ret < 0)
    perror("writing to log file");
}

bool FakeRateLooper::fillErrorInPrediction(TH1F* prediction,
					   TH3F* predictionError,
					   bool addStatisticalError) {
  //
  // calculate error for prediction from predictionError
  //

  for (int predictionBin = 0;
       predictionBin <= prediction->GetNbinsX()+1;
       ++predictionBin ) {
    float err2 = 0;
    for ( int fakeXBin = 0;
	  fakeXBin <= fakeRate().GetNbinsX()+1;
	  ++fakeXBin ) {
      for ( int fakeYBin = 0;
	    fakeYBin <= fakeRate().GetNbinsY()+1;
	    ++fakeYBin ) {
	err2 += predictionError->GetBinContent(predictionBin,fakeXBin,fakeYBin) * 
	  predictionError->GetBinContent(predictionBin,fakeXBin,fakeYBin);
      }
    }
    float err = 0;
    if ( addStatisticalError ) {
      err = sqrt( prediction->GetBinError(predictionBin) *
		  prediction->GetBinError(predictionBin) +
		  err2 );
    } else {
      err = sqrt(err2);
    }

    prediction->SetBinError(predictionBin,err);
  }
  return true;
}

bool FakeRateLooper::extractBins(TAxis *axis, unsigned int nBins, float* bins) {
  //
  // extract float array of bin borders with dimension nbins+1
  //

  for ( unsigned int bin = 0;
	bin < nBins+1;
	++bin ) {
    bins[bin] = axis->GetBinLowEdge(bin+1);
  }
  return true;
}

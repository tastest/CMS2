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
  memset(cands_passing_syst_hi, 0, sizeof(cands_passing_syst_hi));
  memset(cands_passing_syst_lo, 0, sizeof(cands_passing_syst_lo));
  memset(cands_passing_event_weight_only_	, 0, sizeof(cands_passing_       ));
  memset(cands_passing_event_weight_only_w2_	, 0, sizeof(cands_passing_w2_    ));
}

void FakeRateLooper::BookHistos 	()
{
  Looper::BookHistos();

  gDirectory = histo_directory;

  fake_syst = new TH2F("fake_syst", "fake syst uncertainty;#eta;pt", 
		       fakeRate().GetNbinsX(), fakeRate().GetXaxis()->GetXbins()->GetArray(),
		       fakeRate().GetNbinsY(), fakeRate().GetYaxis()->GetXbins()->GetArray());

  const unsigned int jetNBins = 5;
  float jetBins[jetNBins+1] = {0.,1.,2.,3.,4.,5.};

  const unsigned int ptNBins = 16;
  float ptBins[ptNBins+1];
  for ( unsigned int ptBin = 0;
	ptBin <= ptNBins;
	++ptBin) {
    ptBins[ptBin] = float(ptBin)*160./16.;
  }

  const unsigned int metNBins = 100;
  float metBins[metNBins+1];
  for ( unsigned int metBin = 0;
	metBin <= metNBins;
	++metBin) {
    metBins[metBin] = float(metBin)*200./100.;
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

  for (unsigned int bucket = 0;
       bucket < 4;
       ++bucket ) {
    hnJet3D_[bucket] = book3DVarHist(Form("%s_%s_%s",sample_.name.c_str(),"nJet3D",dilepton_hypo_names[bucket]),
				     Form("%s_%s_%s",sample_.name.c_str(),"nJet3D",dilepton_hypo_names[bucket]),
				     jetNBins,jetBins,
				     fakeXNBins,fakeXBins,
				     fakeYNBins,fakeYBins,
				     "n_{jet}","#eta", "p_{T} [GeV]", sample_.histo_color);
    helPt3D_[bucket] = book3DVarHist(Form("%s_%s_%s",sample_.name.c_str(),"elPt3D",dilepton_hypo_names[bucket]),
				     Form("%s_%s_%s",sample_.name.c_str(),"elPt3D",dilepton_hypo_names[bucket]),
				     ptNBins,ptBins,
				     fakeXNBins,fakeXBins,
				     fakeYNBins,fakeYBins,
				     "p_{T}^{e} [GeV]","#eta", "p_{T} [GeV]",sample_.histo_color);
    helEta3D_[bucket] = book3DVarHist(Form("%s_%s_%s",sample_.name.c_str(),"elEta3D",dilepton_hypo_names[bucket]),
				      Form("%s_%s_%s",sample_.name.c_str(),"elEta3D",dilepton_hypo_names[bucket]),
				      etaNBins,etaBins,
				      fakeXNBins,fakeXBins,
				      fakeYNBins,fakeYBins,
				      "#eta^{e} [GeV]","#eta", "p_{T} [GeV]",sample_.histo_color);
    hmet3D_[bucket] = book3DVarHist(Form("%s_%s_%s",sample_.name.c_str(),"met3D",dilepton_hypo_names[bucket]),
				    Form("%s_%s_%s",sample_.name.c_str(),"met3D",dilepton_hypo_names[bucket]),
				    metNBins,metBins,
				    fakeXNBins,fakeXBins,
				    fakeYNBins,fakeYBins,
				    "MET [GeV]","#eta", "p_{T} [GeV]",sample_.histo_color);
  }

}

cuts_t FakeRateLooper::DilepSelect (int i_hyp)
{
  // this doesn't do anything special
  return Looper::DilepSelect(i_hyp);
}

void FakeRateLooper::FillDilepHistos (int i_hyp)
{
  const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);
  cuts_t cuts_passed = DilepSelect(i_hyp);
     
  if ((cuts_passed & cuts_) == cuts_) {

    const double weight_hi = Weight(i_hyp, 1);
    const double weight_lo = Weight(i_hyp, -1);
    cands_passing_syst_hi[myType] += weight_hi;
    cands_passing_syst_lo[myType] += weight_lo;
    cands_passing_syst_hi[DILEPTON_ALL] += weight_hi;
    cands_passing_syst_lo[DILEPTON_ALL] += weight_lo;
    const double weight = Looper::Weight(i_hyp);
    cands_passing_event_weight_only_[myType] += weight;
    cands_passing_event_weight_only_w2_[myType] += weight * weight;
    cands_passing_event_weight_only_[DILEPTON_ALL] += weight;
    cands_passing_event_weight_only_w2_[DILEPTON_ALL] += weight * weight;

    // this doesn't work for ee, since it assumes only one of
    // the hyp leptons is an electron
    if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
      const double err = elFakeProb(cms2.hyp_lt_index()[i_hyp], 1) - 
	elFakeProb(cms2.hyp_lt_index()[i_hyp], 0);
      const double eta = cms2.els_p4()[cms2.hyp_lt_index()[i_hyp]].eta();
      const double pt = cms2.els_p4()[cms2.hyp_lt_index()[i_hyp]].pt();
      fake_syst->Fill(eta, pt, weight * err);
      hnJet3D_[myType]->Fill(cms2.hyp_njets()[i_hyp],eta,pt, weight * err);
      helPt3D_[myType]->Fill(cms2.hyp_lt_p4()[i_hyp].pt(),eta,pt, weight * err);
      helEta3D_[myType]->Fill(cms2.hyp_lt_p4()[i_hyp].eta(),eta,pt, weight * err);
      hmet3D_[myType]->Fill(cms2.hyp_met()[i_hyp],eta,pt, weight * err);
    } else if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
      const double err = elFakeProb(cms2.hyp_ll_index()[i_hyp], 1) - 
	elFakeProb(cms2.hyp_ll_index()[i_hyp], 0);
      const double eta = cms2.els_p4()[cms2.hyp_ll_index()[i_hyp]].eta();
      const double pt = cms2.els_p4()[cms2.hyp_ll_index()[i_hyp]].pt();
      fake_syst->Fill(eta, pt, weight * err);
      hnJet3D_[myType]->Fill(cms2.hyp_njets()[i_hyp],eta,pt, weight * err);
      helPt3D_[myType]->Fill(cms2.hyp_ll_p4()[i_hyp].pt(),eta,pt, weight * err);
      helEta3D_[myType]->Fill(cms2.hyp_ll_p4()[i_hyp].eta(),eta,pt, weight * err);
      hmet3D_[myType]->Fill(cms2.hyp_met()[i_hyp],eta,pt, weight * err);
    }

  }

  Looper::FillDilepHistos(i_hyp);
}

double FakeRateLooper::Weight (int i_hyp)
{ 
  return Weight(i_hyp, 0);
}

double FakeRateLooper::Weight (int i_hyp, int n_sig_syst)
{
  double weight = Looper::Weight(i_hyp);
  double fr = 0;
  // this doesn't work for ee, since it assumes only one of
  // the hyp leptons is an electron
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
    fr = elFakeProb(cms2.hyp_lt_index()[i_hyp], n_sig_syst);
  } else if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
    fr = elFakeProb(cms2.hyp_ll_index()[i_hyp], n_sig_syst);
  }
  return weight * fr / (1 - fr); 
}

double FakeRateLooper::FakeSyst (enum DileptonHypType i) const
{
  switch (i) {
  case DILEPTON_EE: case DILEPTON_EMU:
    return 0;
  default:
    break;
  }
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
  for (unsigned int bucket = 0;
       bucket < 4;
       ++bucket ) {
    fillErrorInPrediction(hnJet_[bucket],hnJet3D_[bucket]);
    fillErrorInPrediction(helPt_[bucket],helPt3D_[bucket]);
    fillErrorInPrediction(helEta_[bucket],helEta3D_[bucket]);
    fillErrorInPrediction(hmet_[bucket],hmet3D_[bucket]);
  }

  ostringstream stream;
  
  stream << endl << "=========" << endl;
  stream << "Sample: " << SampleName().c_str() << endl;
  stream << "=========" << endl;
  stream << "Total candidate count ee: " << CandsCount(DILEPTON_EE) 
	 << " em: " << CandsCount(DILEPTON_EMU) 
	 << " mm: " << CandsCount(DILEPTON_MUMU)
	 << " all: " << CandsCount(DILEPTON_ALL) << endl;
  stream << "Total weighted events (event weight only): ee: " << fixed << setprecision(1) << CandsPassingEventWeightOnly(DILEPTON_EE) << "+-" << RMSEventWeightOnly(DILEPTON_EE)
	 << " em: " << CandsPassingEventWeightOnly(DILEPTON_EMU) << "+-" << RMSEventWeightOnly(DILEPTON_EMU)
	 << " mm: " << CandsPassingEventWeightOnly(DILEPTON_MUMU) << "+-" << RMSEventWeightOnly(DILEPTON_MUMU)
	 << " all: " << CandsPassingEventWeightOnly(DILEPTON_ALL) << "+-" << RMSEventWeightOnly(DILEPTON_ALL) << endl;
  stream << "Total weighted events (weight*fr/(1-fr)): ee: " << fixed << setprecision(1) << CandsPassing(DILEPTON_EE) << "+-" << RMS(DILEPTON_EE)
	 << " em: " << CandsPassing(DILEPTON_EMU) << "+-" << RMS(DILEPTON_EMU)
	 << " mm: " << CandsPassing(DILEPTON_MUMU) << "+-" << RMS(DILEPTON_MUMU)
	 << " all: " << CandsPassing(DILEPTON_ALL) << "+-" << RMS(DILEPTON_ALL) << endl;
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


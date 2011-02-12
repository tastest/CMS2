//
// Studies:
// * top tagging efficiency estimation
//   * Monte Carlo truth (ttbar+tw)
//   * 1-jet
//     * Monte Carlo 
//     * Data
//   * 1b-jet
//     * Monte Carlo 
//     * Data
// * top background estimation in data and Monte Carlo
//

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSystem.h"
#include <iostream>
#include <iomanip>
#include <string>
#include "wwtypes.cc"

struct Efficiency{
  double value;
  double error;
  Efficiency():value(0),error(0){}
};

struct SampleHist{
  TH2F* ee;
  TH2F* em;
  TH2F* mm;
  TH2F* all;
  SampleHist(TFile* ftt, TString pattern){
    ee = dynamic_cast<TH2F*>(ftt->Get(pattern+"_ee"));
    em = dynamic_cast<TH2F*>(ftt->Get(pattern+"_em"));
    mm = dynamic_cast<TH2F*>(ftt->Get(pattern+"_mm"));
    all = dynamic_cast<TH2F*>(ftt->Get(pattern+"_all"));
  }
  bool valid() {return ee && em && mm && all;}
};

enum EvtType { Nominal, OneJet, OneBJet };

Efficiency getTaggingEfficiency(std::vector<SampleHist> samples, EvtType type, bool binomialError=false)
{
  Efficiency eff;
  unsigned int bin(0);
  switch (type){
  case Nominal:
    bin = 1;
    break;
  case OneJet:
    bin = 2;
    break;
  case OneBJet:
    bin = 3;
    break;
  }
  double A(0);
  double B(0);
  double err2A(0);
  double err2B(0);
  for ( std::vector<SampleHist>::const_iterator sample = samples.begin(); 
	sample != samples.end(); ++sample){
    A += sample->all->Integral(2,4,bin,bin);
    B += sample->all->Integral(1,4,bin,bin);
    double x = pow(sample->all->GetBinError(2,bin),2) + pow(sample->all->GetBinError(3,bin),2) +
      pow(sample->all->GetBinError(4,bin),2);
    err2A += x;
    err2B += x + pow(sample->all->GetBinError(1,bin),2);
  }
  if ( B<=0 ) return eff;
  eff.value = A/B;
  if ( binomialError ){
    eff.error = sqrt(eff.value*(1-eff.value)/B);
  } else {
    eff.error = eff.value>0 ? eff.value * sqrt(err2A/pow(A,2)+err2B/pow(B,2)): 0;
  }
  return eff;
}

Efficiency getTaggingEfficiency(SampleHist sample, EvtType type, bool binomialError=false)
{
  std::vector<SampleHist> samples;
  samples.push_back(sample);
  return getTaggingEfficiency(samples,type,binomialError);
}

Efficiency getDoubleTaggableEfficiencyEstimate(Efficiency singleTaggableEfficiency){
  Efficiency eff;
  eff.value = 1-pow(1-singleTaggableEfficiency.value,2);
  eff.error = 2*(1-singleTaggableEfficiency.value)*singleTaggableEfficiency.error;
  return eff;
}

void estimateTopBkgOnZ(const char* file = "processed_data.root")
{
  using namespace std;
  
  //hist file:
  TFile *ftt = TFile::Open(file);
  assert(ftt);
  
  SampleHist DYee(ftt,"dyee_htoptagz"); 
  SampleHist DYmm(ftt,"dymm_htoptagz"); 
  SampleHist DYtt(ftt,"dytt_htoptagz"); 
  SampleHist tt(ftt,"ttbar_htoptagz");
  SampleHist tw(ftt,"tw_htoptagz");
  SampleHist wjets(ftt,"wjets_htoptagz");
  SampleHist wz(ftt,"wz_htoptagz");
  SampleHist zz(ftt,"zz_htoptagz");
  SampleHist ww(ftt,"ww_htoptagz");
  SampleHist data(ftt,"data_htoptagz");

  if ( DYee.valid() && DYmm.valid() &&
       ww.valid() && data.valid() ) {
    Efficiency ww_nom = getTaggingEfficiency(ww,Nominal);

    std::vector<SampleHist> samples;
    samples.push_back(DYee);
    samples.push_back(DYmm);
    Efficiency dy_nom = getTaggingEfficiency(samples,Nominal);
    Efficiency data_nom = getTaggingEfficiency(data,Nominal);

    printf("\nZ-sample for miss-tag:\n");
    printf("WW MC (nominal): %0.2f%% +/- %0.2f%%\n", ww_nom.value*100, ww_nom.error*100) ;
    printf("DY MC (nominal): %0.2f%% +/- %0.2f%%\n", dy_nom.value*100, dy_nom.error*100) ;
    printf("Data in Z-window without MET requirement (nominal): %0.2f%% +/- %0.2f%%\n", data_nom.value*100, data_nom.error*100) ;
  }
}

void estimateTopBkg(const char* file = "processed_data.root")
{
  using namespace std;
  
  //hist file:
  TFile *ftt = TFile::Open(file);
  assert(ftt);
  
  SampleHist DYee(ftt,"dyee_htoptag"); 
  SampleHist DYmm(ftt,"dymm_htoptag"); 
  SampleHist DYtt(ftt,"dytt_htoptag"); 
  SampleHist tt(ftt,"ttbar_htoptag");
  SampleHist tw(ftt,"tw_htoptag");
  SampleHist wjets(ftt,"wjets_htoptag");
  SampleHist wz(ftt,"wz_htoptag");
  SampleHist zz(ftt,"zz_htoptag");
  SampleHist ww(ftt,"ww_htoptag");
  SampleHist data(ftt,"data_htoptag");

  if ( tt.valid() && tw.valid() ){
    // Sanity check first
    TH2F* ttbar_toptagvsnjet_all = dynamic_cast<TH2F*>(ftt->Get("ttbar_toptagvsnjet_all"));
    if ( ttbar_toptagvsnjet_all ){
      double a = ttbar_toptagvsnjet_all->Integral(3,10,1,1);
      double b = ttbar_toptagvsnjet_all->Integral(3,10,2,2);
      printf("ttbar MC (2 or more jets) top tagging efficiency: %0.1f%%\n", a+b>0?100*b/(a+b):0);
    }
    TH2F* data_toptagvsnjet_all = dynamic_cast<TH2F*>(ftt->Get("data_toptagvsnjet_all"));
    if ( ttbar_toptagvsnjet_all ){
      double a = data_toptagvsnjet_all->Integral(3,10,1,1);
      double b = data_toptagvsnjet_all->Integral(3,10,2,2);
      printf("ttbar data (2 or more jets) top tagging efficiency: %0.1f%%\n\n", a+b>0?100*b/(a+b):0);
    }
    Efficiency tt_nom = getTaggingEfficiency(tt,Nominal);
    Efficiency tt_1jet = getTaggingEfficiency(tt,OneJet);
    Efficiency tt_1bjet = getTaggingEfficiency(tt,OneBJet);
    Efficiency tt_1jet_est = getDoubleTaggableEfficiencyEstimate(tt_1jet);
    Efficiency tt_1bjet_est = getDoubleTaggableEfficiencyEstimate(tt_1bjet);

    Efficiency tw_nom = getTaggingEfficiency(tw,Nominal);
    Efficiency tw_1jet = getTaggingEfficiency(tw,OneJet);
    Efficiency tw_1bjet = getTaggingEfficiency(tw,OneBJet);
    Efficiency tw_1jet_est = getDoubleTaggableEfficiencyEstimate(tw_1jet);
    Efficiency tw_1bjet_est = getDoubleTaggableEfficiencyEstimate(tw_1bjet);
    
    std::vector<SampleHist> samples;
    samples.push_back(tt);
    samples.push_back(tw);
    Efficiency top_nom = getTaggingEfficiency(samples,Nominal);
    Efficiency top_1jet = getTaggingEfficiency(samples,OneJet);
    Efficiency top_1bjet = getTaggingEfficiency(samples,OneBJet);
    Efficiency top_1jet_est = getDoubleTaggableEfficiencyEstimate(top_1jet);
    Efficiency top_1bjet_est = getDoubleTaggableEfficiencyEstimate(top_1bjet);

    printf("Monte Carlo truth based tagging efficiency:\n");
    printf("ttbar (nominal): %0.2f%% +/- %0.2f%%\n",tt_nom.value*100, tt_nom.error*100) ;
    printf("\tttbar (1jet): %0.2f%% +/- %0.2f%% \test(2jet): %0.2f%% +/- %0.2f%%\n",
	   100*tt_1jet.value, 100*tt_1jet.error, 100*tt_1jet_est.value, 100*tt_1jet_est.error);
    printf("\tttbar (1bjet): %0.2f%% +/- %0.2f%% \test(2jet): %0.2f%% +/- %0.2f%%\n",
	   100*tt_1bjet.value, 100*tt_1bjet.error, 100*tt_1bjet_est.value, 100*tt_1bjet_est.error);
    
    printf("tw (nominal): %0.2f%% +/- %0.2f%%\n",tw_nom.value*100, tw_nom.error*100) ;
    printf("\ttw (1jet): %0.2f%% +/- %0.2f%% \test(2jet): %0.2f%% +/- %0.2f%%\n",
	   100*tw_1jet.value, 100*tw_1jet.error, 100*tw_1jet_est.value, 100*tw_1jet_est.error);
    printf("\ttw (1bjet): %0.2f%% +/- %0.2f%% \test(2jet): %0.2f%% +/- %0.2f%%\n",
	   100*tw_1bjet.value, 100*tw_1bjet.error, 100*tw_1bjet_est.value, 100*tw_1bjet_est.error);
    
    printf("top (nominal): %0.2f%% +/- %0.2f%%\n",top_nom.value*100, top_nom.error*100) ;
    printf("\ttop (1jet): %0.2f%% +/- %0.2f%% \test(2jet): %0.2f%% +/- %0.2f%%\n",
	   100*top_1jet.value, 100*top_1jet.error, 100*top_1jet_est.value, 100*top_1jet_est.error);
    printf("\ttop (1bjet): %0.2f%% +/- %0.2f%% \test(2jet): %0.2f%% +/- %0.2f%%\n",
	   100*top_1bjet.value, 100*top_1bjet.error, 100*top_1bjet_est.value, 100*top_1bjet_est.error);

    if ( DYee.valid() && DYmm.valid() && DYtt.valid() &&
	 wjets.valid() && wz.valid() && zz.valid() && ww.valid() ){
      samples.push_back(DYee);
      samples.push_back(DYmm);
      samples.push_back(DYtt);
      samples.push_back(wjets);
      samples.push_back(wz);
      samples.push_back(zz);
      samples.push_back(ww);
      Efficiency all_1jet = getTaggingEfficiency(samples,OneJet);
      Efficiency all_1bjet = getTaggingEfficiency(samples,OneBJet);
      Efficiency all_1jet_est = getDoubleTaggableEfficiencyEstimate(all_1jet);
      Efficiency all_1bjet_est = getDoubleTaggableEfficiencyEstimate(all_1bjet);

      printf("All (1jet): %0.2f%% +/- %0.2f%% \test(2jet): %0.2f%% +/- %0.2f%%\n",
	     100*all_1jet.value, 100*all_1jet.error, 100*all_1jet_est.value, 100*all_1jet_est.error);
      printf("All (1bjet): %0.2f%% +/- %0.2f%% \test(2jet): %0.2f%% +/- %0.2f%%\n",
	     100*all_1bjet.value, 100*all_1bjet.error, 100*all_1bjet_est.value, 100*all_1bjet_est.error);
    }
    if ( data.valid() ){
      Efficiency data_1jet = getTaggingEfficiency(data,OneJet,true);
      Efficiency data_1jet_2 = data_1jet; 
      data_1jet_2.error = sqrt(top_1jet.value*(1-top_1jet.value)/data.all->Integral(1,4,2,2));
      Efficiency data_1bjet = getTaggingEfficiency(data,OneBJet,true);
      Efficiency data_1bjet_2 = data_1bjet; 
      data_1bjet_2.error = sqrt(top_1bjet.value*(1-top_1bjet.value)/data.all->Integral(1,4,2,2));
      Efficiency data_1jet_est = getDoubleTaggableEfficiencyEstimate(data_1jet);
      Efficiency data_1jet_est_2 = getDoubleTaggableEfficiencyEstimate(data_1jet_2);
      Efficiency data_1bjet_est = getDoubleTaggableEfficiencyEstimate(data_1bjet);
      Efficiency data_1bjet_est_2 = getDoubleTaggableEfficiencyEstimate(data_1bjet_2);
      printf("Data (1jet): %0.2f%% +/- %0.2f%% (+/- %0.2f%%) \test(2jet): %0.2f%% +/- %0.2f%% (+/- %0.2f%%)\n",
	     100*data_1jet.value, 100*data_1jet.error, 100*data_1jet_2.error,
	     100*data_1jet_est.value, 100*data_1jet_est.error, 100*data_1jet_est_2.error);
      printf("Data (1bjet): %0.2f%% +/- %0.2f%% (+/- %0.2f%%) \test(2jet): %0.2f%% +/- %0.2f%% (+/- %0.2f%%)\n",
	     100*data_1bjet.value, 100*data_1bjet.error, 100*data_1bjet_2.error,
	     100*data_1bjet_est.value, 100*data_1bjet_est.error, 100*data_1bjet_est_2.error);
    }
  } else {
    cout << "top samples are not available" << endl;
  }
  estimateTopBkgOnZ(file);
}


#include "cuts.h"
#include "gather.h"
#include "goodrun.h"

#include "TCanvas.h"
#include "TChain.h"
#include "TH1F.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include <iostream>
#include <vector>

TH1F * slideIntegrated(TH1F* integrateThis) {
  //  std::cout << integrateThis->GetName() << std::endl;
  TString name = integrateThis->GetName();
  int NbinsX = integrateThis->GetNbinsX();
  TH1F* integrated = ((TH1F*)integrateThis->Clone("integrated"));
  integrated->Reset();

  // scale histo to account for the number of simulated muons
  // 2*10^34, pt 1-125  eta 1.2-1.9 
  //  integrateThis->Scale(2*124);
  // integrate rate above a certain pt value (bin)
  for(int i = 1; i<NbinsX; ++i) {
    float integral = integrateThis->Integral(i,NbinsX);
    if(integral!=0.) integrated->SetBinContent(i,integral);
    if(integral!=0.) integrated->SetBinError(i,0);
  }

  TH1F * result = ((TH1F*)integrated->Clone("name"));  
  //  result->Smooth(80);

  return result;
}

void  DrawAll( TString drawThis,TCut thisBaseSel,float f_intlumifb,float f_kfactor, UInt_t nBinsX, float lowBinX, float highBinX, Bool_t  integrated, TCanvas * c1,
               TChain * cdata_dilep      ,
               TChain * cqcd15_dilep     ,
               TChain * cqcd30_dilep     ,
               TChain * cqcd80_dilep     ,
               TChain * cqcd170_dilep    ,
               TChain * cttbar_dilep     ,
               TChain * cttbarjets_dilep ,
               TChain * csingletop_dilep ,
               TChain * cdyee_dilep      ,
               TChain * cdymumu_dilep    ,
               TChain * cvvjets_dilep    ,
               TChain * cwjets_dilep     ,
               TChain * czjets_dilep     ,
               TChain * czee_dilep       ,
               TChain * czmumu_dilep     ,
               TChain * cztautau_dilep  
               ) {

  std::vector<TH1F*> hmcs;
  TH1F * buffer;
      
  TLegend* leg = new TLegend(0.65,0.622881,0.997126,0.98,NULL,"brNDC");
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(10);
  leg->SetBorderSize(1);
      
      
  if(42 != 42) {
    buffer = Plot("ttbar"     , cttbar_dilep     ,drawThis,thisBaseSel,"",f_intlumifb,f_kfactor,nBinsX,lowBinX,highBinX);             hmcs.push_back((TH1F*)buffer->Clone("hmc_ttbar"));
  }
  buffer = Plot("ttbarjets" , cttbarjets_dilep ,drawThis,thisBaseSel,"",f_intlumifb,f_kfactor,nBinsX,lowBinX,highBinX);             hmcs.push_back((TH1F*)buffer->Clone("hmc_ttbarjets"));
  buffer = Plot("singletop" , csingletop_dilep ,drawThis,thisBaseSel,"",f_intlumifb,f_kfactor,nBinsX,lowBinX,highBinX);             hmcs.push_back((TH1F*)buffer->Clone("hmc_singletop"));
      
  buffer = Plot("vvjets"    , cvvjets_dilep    ,drawThis,thisBaseSel,"",f_intlumifb,f_kfactor,nBinsX,lowBinX,highBinX);             hmcs.push_back((TH1F*)buffer->Clone("hmc_vvjets"));
  buffer = Plot("wjets"     , cwjets_dilep     ,drawThis,thisBaseSel,"",f_intlumifb,f_kfactor,nBinsX,lowBinX,highBinX);             hmcs.push_back((TH1F*)buffer->Clone("hmc_wjets"));
  // m_ll > 50 GeV
  buffer = Plot("zjets"     , czjets_dilep     ,drawThis,thisBaseSel,"",f_intlumifb,f_kfactor,nBinsX,lowBinX,highBinX);             hmcs.push_back((TH1F*)buffer->Clone("hmc_zjets"));
      
  // m_ll < 50 GeV
  buffer = Plot("zee"       , czee_dilep       ,drawThis,thisBaseSel,"mass < 50.",f_intlumifb,f_kfactor,nBinsX,lowBinX,highBinX);   hmcs.push_back((TH1F*)buffer->Clone("hmc_zee"));
  buffer = Plot("zmumu"     , czmumu_dilep     ,drawThis,thisBaseSel,"mass < 50.",f_intlumifb,f_kfactor,nBinsX,lowBinX,highBinX);   hmcs.push_back((TH1F*)buffer->Clone("hmc_zmumu"));
  buffer = Plot("ztautau"   , cztautau_dilep   ,drawThis,thisBaseSel,"mass < 50.",f_intlumifb,f_kfactor,nBinsX,lowBinX,highBinX);   hmcs.push_back((TH1F*)buffer->Clone("hmc_ztautau"));
  // m_ll 10 - 20 GeV                                                                                        
  buffer = Plot("dyee"      , cdyee_dilep      ,drawThis,thisBaseSel,"",f_intlumifb,f_kfactor,nBinsX,lowBinX,highBinX);             hmcs.push_back((TH1F*)buffer->Clone("hmc_dyee"));
  buffer = Plot("dymumu"    , cdymumu_dilep    ,drawThis,thisBaseSel,"",f_intlumifb,f_kfactor,nBinsX,lowBinX,highBinX);             hmcs.push_back((TH1F*)buffer->Clone("hmc_dymumu"));
      
  buffer = Plot("qcd15",    cqcd15_dilep     ,drawThis,thisBaseSel,"pthat<30",f_intlumifb,f_kfactor,nBinsX,lowBinX,highBinX);       hmcs.push_back((TH1F*)buffer->Clone("hmc_qcd15"));
  buffer = Plot("qcd30",    cqcd30_dilep     ,drawThis,thisBaseSel,"pthat<80",f_intlumifb,f_kfactor,nBinsX,lowBinX,highBinX);       hmcs.push_back((TH1F*)buffer->Clone("hmc_qcd30"));
  buffer = Plot("qcd80",    cqcd80_dilep     ,drawThis,thisBaseSel,"pthat<170",f_intlumifb,f_kfactor,nBinsX,lowBinX,highBinX);      hmcs.push_back((TH1F*)buffer->Clone("hmc_qcd80"));
  buffer = Plot("qcd170",   cqcd170_dilep    ,drawThis,thisBaseSel,"",f_intlumifb,f_kfactor,nBinsX,lowBinX,highBinX);               hmcs.push_back((TH1F*)buffer->Clone("hmc_qcd170"));                     
      
  // TH1F* Plot(const char *pfx, TChain *chain, const char *field, TCut sel, TCut presel,
  //                   unsigned int nbins, float xlo, float xhi);
  // same as above but as this is data there is no need to specify intlumifb or kfactor
  TH1F* hdata = Plot("data", cdata_dilep, drawThis,thisBaseSel,"",nBinsX,lowBinX,highBinX);
      
  // there is another Plot implementation for data as well:
  // TH1F* Plot(const char *pfx, TChain *chain, const char *field, TCut sel, TCut presel, float kfactor,
  //                   unsigned int nbins, float xlo, float xhi);
  // which is the same as above but allows specification of a kfactor, which I can imagine using
      
  // mini Hstack function: (should be replaced with something more adequate)
  for(uint mchisto = 0; mchisto < hmcs.size(); ++mchisto) {
    //      std::cout<<"We are at histo "<<mchisto<<std::endl;
    for(uint mchisto2 = mchisto+1; mchisto2 < hmcs.size(); ++mchisto2) {
      //        std::cout<<"Adding histos "<<mchisto2<<std::endl;
      (hmcs.at(mchisto))->Add(hmcs.at(mchisto2),1);
    }
  }
      
  if(integrated)    hmcs[0] = slideIntegrated( hmcs[0] );
  if(integrated)    hdata = slideIntegrated( hdata );
      
  float ymax = hdata->GetMaximum() > hmcs[0]->GetMaximum() ? hdata->GetMaximum() : hmcs[0]->GetMaximum();
  hdata->SetMarkerStyle(20);
  leg->AddEntry(hdata, hdata->GetTitle(),"lep");
      
  hmcs[0]->SetFillColor(kOrange);
  hmcs[0]->SetLineColor(kOrange);
  hmcs[0]->SetMaximum(ymax*1.1);
  // Draw the initial MC histo
  hmcs[0]->Draw();
  leg->AddEntry(hmcs.at(0), hmcs.at(0)->GetTitle(),"f");
      
  for(uint mchisto = 1; mchisto < hmcs.size(); ++mchisto) { // need 1 here, as we only want to plot 1 and more
    hmcs.at(mchisto)->SetFillColor(mchisto+1);
    hmcs.at(mchisto)->SetLineColor(mchisto+1);
        
    if(integrated)      hmcs.at(mchisto) = slideIntegrated( hmcs.at(mchisto) );
        
    hmcs.at(mchisto)->Draw("same");
        
    leg->AddEntry(hmcs.at(mchisto), hmcs.at(mchisto)->GetTitle(),"f");
  }
      
  hdata->Draw("pesames");
      
  leg->Draw();
      
  c1->RedrawAxis();
}

void gather_several_doAll()
{
  gROOT->SetStyle("Plain");
  gStyle->SetHistMinimumZero();
  //  gStyle->SetOptStat("nemoui");
  gStyle->SetOptStat(0);
  
   // this sets the json file obviously; it's
    // easiest to just  soft link to  the real
    // file as json.txt
    std::cout << "Using json.txt for goodruns\n";
    set_goodrun_file_json("json.txt");

    // calculate  integrated luminosity  in /fb
    // note that I input /nb and so get out /nb
    // then scale to /fb which  is what  I need
    // 1e-3*GetIntLumi(/pb)  is the  same thing
    // doesn't matter so long as it corresponds
    // to the lumi of the json file and that it
    // is correctly scaled to /fb afterward
    //    float f_intlumifb = 1e-6*GetIntLumi(2790);
    float f_intlumifb = 1e-6*GetIntLumi(10600); // random guess to make the MC get to the data...
    std::cout << "Integrated luminosity: " << f_intlumifb << "/fb\n";


    // data babies
    TChain *cdata_emu = new TChain("tree");
    cdata_emu->Add("/tas05/disk00/jribnik/hunt/emu_baby/*.root");
    TChain *cdata_dilep = new TChain("tree");
    cdata_dilep->Add("/tas05/disk00/jribnik/hunt/dilep_baby/*.root");
    TChain *cdata_trilep = new TChain("tree");
    cdata_trilep->Add("/tas05/disk00/jribnik/hunt/trilep_baby/*.root");

    // -----------------
    // mc babies
    // -----------------
    //
    // singlelep babies
    //
    TChain *czjets_emu = new TChain("tree");
    czjets_emu->Add("/tas05/disk00/jribnik/huntmc/ZJets-madgraph_Spring10-START3X_V26_S09-v1/emu_baby/*.root");
    //
    //dilep babies
    //
    TChain *cqcd15_dilep     = new TChain("tree");
    TChain *cqcd30_dilep     = new TChain("tree");
    TChain *cqcd80_dilep     = new TChain("tree");
    TChain *cqcd170_dilep    = new TChain("tree");

    TChain *cttbar_dilep     = new TChain("tree");
    TChain *cttbarjets_dilep = new TChain("tree");
    TChain *csingletop_dilep = new TChain("tree");

    TChain *cdyee_dilep      = new TChain("tree");
    TChain *cdymumu_dilep    = new TChain("tree");

    TChain *cvvjets_dilep    = new TChain("tree");
    TChain *cwjets_dilep     = new TChain("tree");
    TChain *czjets_dilep     = new TChain("tree");

    TChain *czee_dilep       = new TChain("tree");
    TChain *czmumu_dilep     = new TChain("tree");
    TChain *cztautau_dilep   = new TChain("tree");

    cqcd15_dilep     ->Add("/tas05/disk00/jribnik/huntmc/QCD_Pt15_Spring10-START3X_V26_S09-v1/dilep_baby/*.root");
    cqcd30_dilep     ->Add("/tas05/disk00/jribnik/huntmc/QCD_Pt30_Spring10-START3X_V26_S09-v1/dilep_baby/*.root");
    cqcd80_dilep     ->Add("/tas05/disk00/jribnik/huntmc/QCD_Pt80_Spring10-START3X_V26_S09-v1/dilep_baby/*.root");
    cqcd170_dilep    ->Add("/tas05/disk00/jribnik/huntmc/QCD_Pt170_Spring10-START3X_V26_S09-v1/dilep_baby/*.root");
    
    cttbar_dilep     ->Add("/tas05/disk00/jribnik/huntmc/TTbar_Spring10-START3X_V26_S09-v1/dilep_baby/*.root");
    cttbarjets_dilep ->Add("/tas05/disk00/jribnik/huntmc/TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root");
    csingletop_dilep ->Add("/tas05/disk00/jribnik/huntmc/SingleTop_tWChannel-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root");
    
    cdyee_dilep      ->Add("/tas05/disk00/jribnik/huntmc/DYee_M10to20_Spring10-START3X_V26_S09-v1/dilep_baby/*.root");
    cdymumu_dilep    ->Add("/tas05/disk00/jribnik/huntmc/DYmumu_M10to20_Spring10-START3X_V26_S09-v1/dilep_baby/*.root");
    
    cvvjets_dilep    ->Add("/tas05/disk00/jribnik/huntmc/VVJets-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root");
    cwjets_dilep     ->Add("/tas05/disk00/jribnik/huntmc/WJets-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root");
    czjets_dilep     ->Add("/tas05/disk00/jribnik/huntmc/ZJets-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root");

    czee_dilep       ->Add("/tas05/disk00/jribnik/huntmc/Zee_Spring10-START3X_V26_S09-v1/dilep_baby/*.root");
    czmumu_dilep     ->Add("/tas05/disk00/jribnik/huntmc/Zmumu_Spring10-START3X_V26_S09-v1/dilep_baby/*.root");
    cztautau_dilep   ->Add("/tas05/disk00/jribnik/huntmc/Ztautau_Spring10-START3X_V26_S09-v1/dilep_baby/*.root");    

    // trilep babies
    //
    TChain *czjets_trilep = new TChain("tree");
    czjets_trilep->Add("/tas05/disk00/jribnik/huntmc/ZJets-madgraph_Spring10-START3X_V26_S09-v1/trilep_baby/*.root");

    // make plots
    // here is a breakdown of the Plot method
    // TH1F* Plot(const char *pfx, TChain *chain, const char *field, TCut sel, TCut presel, float intlumifb, float kfactor,
    //                   unsigned int nbins, float xlo, float xhi);
    // the arguments are:
    //     prefix for histogram names
    //     chain of files
    //     the baby-ntuple-field you want to draw
    //     a selection, generally one used by the hunt defined in cuts.h
    //     a pre-selection, this is applied in addition to the selection, it is meant for
    //         things like stitching together MC, i.e. "pthat<30"; I imagine it can  also
    //         be used to add OS,SS,OF,SF cuts, i.e. "abs(eormu1)!=abs(eormu2) w/o having
    //         to create a new TCut in cuts.h
    //     the integrated luminosity in /fb; scale1fb is of course multiplied by this no.
    //     whatever additional scale factor you need
    //     then the standard nbins and histogram range

    // prep
    // TCanvas * ttbarPlots(drawThis,thisBaseSel,nBinsX,60,120);
    // enum ttbar_enum { ttbar = 0, dymm, ... }
    // float ttbar_kfactors[sizeofttbarsampels] 

    Bool_t  integrated  = false;
    TString drawThis    = "mass";
    TCut    thisBaseSel = zmet30_dilep;
    UInt_t  nBinsX      = 100;
    Float_t lowBinX     = 0.;
    Float_t highBinX    = 120.;

    TCanvas * c1 = new TCanvas("c1","c1");

    DrawAll( "mass", thisBaseSel, f_intlumifb, 1., 100, 70, 170, integrated, c1,
             cdata_dilep      ,
             cqcd15_dilep     ,
             cqcd30_dilep     ,
             cqcd80_dilep     ,
             cqcd170_dilep    ,
             cttbar_dilep     ,
             cttbarjets_dilep ,
             csingletop_dilep ,
             cdyee_dilep      ,
             cdymumu_dilep    ,
             cvvjets_dilep    ,
             cwjets_dilep     ,
             czjets_dilep     ,
             czee_dilep       ,
             czmumu_dilep     ,
             cztautau_dilep  
             );
    c1->SaveAs("mass.png");

    TCanvas * c2 = new TCanvas("c2","c2");
    DrawAll( "mt2", thisBaseSel, f_intlumifb, 1., 100, 0, 120, true, c2,
             cdata_dilep      ,
             cqcd15_dilep     ,
             cqcd30_dilep     ,
             cqcd80_dilep     ,
             cqcd170_dilep    ,
             cttbar_dilep     ,
             cttbarjets_dilep ,
             csingletop_dilep ,
             cdyee_dilep      ,
             cdymumu_dilep    ,
             cvvjets_dilep    ,
             cwjets_dilep     ,
             czjets_dilep     ,
             czee_dilep       ,
             czmumu_dilep     ,
             cztautau_dilep  
             );
    c2->SaveAs("mt2_int.png");


    TCanvas * c3 = new TCanvas("c3","c3");
    DrawAll( "tcmet", thisBaseSel, f_intlumifb, 1., 100, 0, 400, true, c3,
             cdata_dilep      ,
             cqcd15_dilep     ,
             cqcd30_dilep     ,
             cqcd80_dilep     ,
             cqcd170_dilep    ,
             cttbar_dilep     ,
             cttbarjets_dilep ,
             csingletop_dilep ,
             cdyee_dilep      ,
             cdymumu_dilep    ,
             cvvjets_dilep    ,
             cwjets_dilep     ,
             czjets_dilep     ,
             czee_dilep       ,
             czmumu_dilep     ,
             cztautau_dilep  
             );
    c3->SaveAs("tcmet_int.png");

    TCanvas * c3a = new TCanvas("c3a","c3a");
    DrawAll( "tcmet", thisBaseSel, f_intlumifb, 1., 100, 0, 400, false, c3a,
             cdata_dilep      ,
             cqcd15_dilep     ,
             cqcd30_dilep     ,
             cqcd80_dilep     ,
             cqcd170_dilep    ,
             cttbar_dilep     ,
             cttbarjets_dilep ,
             csingletop_dilep ,
             cdyee_dilep      ,
             cdymumu_dilep    ,
             cvvjets_dilep    ,
             cwjets_dilep     ,
             czjets_dilep     ,
             czee_dilep       ,
             czmumu_dilep     ,
             cztautau_dilep  
             );
    c3a->SaveAs("tcmet.png");

    TCanvas * c4 = new TCanvas("c4","c4");
    DrawAll( "njetsClean", thisBaseSel, f_intlumifb, 1., 15, 0, 15, false, c4,
             cdata_dilep      ,
             cqcd15_dilep     ,
             cqcd30_dilep     ,
             cqcd80_dilep     ,
             cqcd170_dilep    ,
             cttbar_dilep     ,
             cttbarjets_dilep ,
             csingletop_dilep ,
             cdyee_dilep      ,
             cdymumu_dilep    ,
             cvvjets_dilep    ,
             cwjets_dilep     ,
             czjets_dilep     ,
             czee_dilep       ,
             czmumu_dilep     ,
             cztautau_dilep  
             );
    c4->SaveAs("njetsClean.png");

    TCanvas * c5 = new TCanvas("c5","c5");
    DrawAll( "dilpt", thisBaseSel, f_intlumifb, 1., 200, 0, 200, false, c5,
             cdata_dilep      ,
             cqcd15_dilep     ,
             cqcd30_dilep     ,
             cqcd80_dilep     ,
             cqcd170_dilep    ,
             cttbar_dilep     ,
             cttbarjets_dilep ,
             csingletop_dilep ,
             cdyee_dilep      ,
             cdymumu_dilep    ,
             cvvjets_dilep    ,
             cwjets_dilep     ,
             czjets_dilep     ,
             czee_dilep       ,
             czmumu_dilep     ,
             cztautau_dilep  
             );
    c5->SaveAs("dilpt_int.png");

 
}


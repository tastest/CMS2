#include <iostream>
#include "histtools.C"
#include "printHists.C"


//void getYields(TString dataFName="results_data/hist_usePtGt2020_applyTriggers_hypDisamb_usejptJets_usetcMET_requireEcalEls_useOS_vetoHypMassLt12_require2BTag_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root", TString mcFName="results/hist_usePtGt2020_applyTriggers_hypDisamb_usejptJets_usetcMET_requireEcalEls_useOS_vetoHypMassLt12_require2BTag_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root", float scaleMC=1.0, bool combineJetBins = false, const char* formatS = "%6.3f", bool latex = true) {
void getYields(TString FName="results/hist_usePtGt2020_applyTriggers_hypDisamb_usepfMET_usepfJets_requireEcalEls_useOS_vetoHypMassLt12_require2BTag_sortJetCandidatesbyDR_applyLeptonJetInvMassCut170_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root", float scaleMC=1.0, bool combineJetBins = false, const char* formatS = "%6.2f", bool latex = true) {


  

    
    
  hist::loadHist(FName.Data(),0,"*_hnBtagJet_allj_*");
    //    hist::loadHist(mcFName.Data(),0,"*_htcmet_allj_*");
    //hist::loadHist(mcFName.Data(),0,"*hdilMass_allj_*");
    //hist::loadHist(mcFName.Data(),0,"*bTagtk*_allj_*");    
    //hist::scale("*_*", scaleMC);
    // TH1F* ttdil= (TH1F*)gFile->Get("ttdil_hnJets");
    //cout << ttdil->Integral()<<endl;
    //hist::scale("ttdil_*", 1.+(9824. - 10086.77)/9344.25); //scale mc@nlo to data
    hist::scale("ttdil_*", 1.+(9824. - 10063.47)/9323.84); //scale mc@nlo to data, after pT reweighting
    hist::loadHist(FName.Data(),0,"data_hnBtagJet_allj_*");
    //hist::loadHist(dataFName.Data(),0,"data_htcmet_allj_*");
    //hist::loadHist(dataFName.Data(),0,"data_*hdilMass_allj_*");
    //hist::loadHist(dataFName.Data(),0,"data_*bTagtk*_allj_*");
    
    
    std::vector<TString> v_prfxsToCombine;
    //    v_prfxsToCombine.push_back("qcd15");
    //v_prfxsToCombine.push_back("qcd30");
    //hist::combineHists(v_prfxsToCombine,"QCD");
    
    
    //v_prfxsToCombine.clear();
    //v_prfxsToCombine.push_back("DYee");
    //v_prfxsToCombine.push_back("DYmm");
    //hist::combineHists(v_prfxsToCombine, "DYeemm");
    
    
    

    
   
    
    printNJets(latex, formatS,"ttprime", false,true,combineJetBins, false,false); 
    //browseStacks( true, false , "ttdil", dataFName, 4, 7, true, false, 3, false, 0, true);
    //browseStacks( true, false , "DYeemm", dataFName, 4, 27, true, false, 3, false, 0, true);
    hist::deleteHistos();
    gDirectory->Clear();
    

    
  
}



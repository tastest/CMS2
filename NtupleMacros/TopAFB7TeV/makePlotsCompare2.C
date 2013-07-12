#include "TH1F.h"
#include "TString.h"
#include "TCanvas.h"
#include <iostream>
#include "TFile.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLatex.h"
#include <fstream>
#include "tdrStyle.C"
#include "CommonFunctions.C"
#include <vector>

using namespace std;

TH1F* histo1;
TH1F* histo2;
//TH1F* histo3;
TH1F* datahisto;


TH1F* AddUpAllMC(TFile* f_mc, char* histname) {

  std::cout << std::endl;
  f_mc->Print();
  std::cout << std::endl;

  std::vector<TString> mcsamplenames;
  mcsamplenames.push_back("ttdil");
  mcsamplenames.push_back("ttotr");
  mcsamplenames.push_back("wjets");
  mcsamplenames.push_back("DYee");
  mcsamplenames.push_back("DYmm");
  mcsamplenames.push_back("DYtautau");
  mcsamplenames.push_back("VV");
  mcsamplenames.push_back("tw");

  std::cout << Form("%s_%s",mcsamplenames[0].Data(),histname) << " with entries " << ((TH1F*)f_mc->Get(Form("%s_%s",mcsamplenames[0].Data(),histname)))->GetEntries() << std::endl;                                                                                                                                                                                                                                                                                     

  TH1F *clone = (TH1F*)f_mc->Get(Form("%s_%s",mcsamplenames[0].Data(),histname))->Clone();
  std::cout << "adding hist " << clone->GetName() << " with entries " << clone->GetEntries() << std::endl;                                                                                                                                                                                                                                                      
  for (unsigned int i = 1; i < mcsamplenames.size(); ++i) {
    TH1F *hist = (TH1F*)f_mc->Get(Form("%s_%s",mcsamplenames[i].Data(),histname));
    std::cout << "adding hist " << hist->GetName() << " with entries " << hist->GetEntries() << std::endl;                                                                                                                                                                                                                                                
    clone->Add(hist);
  }

  return clone;
}



//void CompareHists(TString histname="hlepAzimAsym2", int rebin = 5, int drawnorm =1, TString FName1 = "/nfs-6/userdata/linacre/AFB_results/no_pT_reweighting/TopAFB_default/results/hist_usePtGt2020_hypDisamb_usepfMET_usepfJets_useOS_vetoHypMassLt12_requireBTag_sortJetCandidatesbyPt_generalLeptonVeto_createBabyNtuples_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root", TString FName2 = "/nfs-6/userdata/linacre/AFB_results/default_madgraph_pT_reweighting/TopAFB_default_madgraph_pT_reweighting/results/hist_usePtGt2020_hypDisamb_usepfMET_usepfJets_useOS_vetoHypMassLt12_requireBTag_sortJetCandidatesbyPt_generalLeptonVeto_createBabyNtuples_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root"){
 void CompareHists(TString histname="httMass", int rebin = 5, int drawnorm =1, TString FName1 = "/nfs-6/userdata/linacre/AFB_results/no_pT_reweighting/TopAFB_default/results/hist_usePtGt2020_hypDisamb_usepfMET_usepfJets_useOS_vetoHypMassLt12_requireBTag_sortJetCandidatesbyPt_generalLeptonVeto_createBabyNtuples_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root", TString FName2 = "/nfs-6/userdata/linacre/AFB_results/default_madgraph_pT_reweighting/TopAFB_default_madgraph_pT_reweighting/results/hist_usePtGt2020_hypDisamb_usepfMET_usepfJets_useOS_vetoHypMassLt12_requireBTag_sortJetCandidatesbyPt_generalLeptonVeto_createBabyNtuples_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root"){
//void CompareHists(TString histname="hlepPt", int rebin = 5, int drawnorm =1, TString FName1 = "/nfs-6/userdata/linacre/AFB_results/no_pT_reweighting/TopAFB_default/results/hist_usePtGt2020_hypDisamb_usepfMET_usepfJets_useOS_vetoHypMassLt12_requireBTag_sortJetCandidatesbyPt_generalLeptonVeto_createBabyNtuples_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root", TString FName2 = "/nfs-6/userdata/linacre/AFB_results/default_madgraph_pT_reweighting/TopAFB_default_madgraph_pT_reweighting/results/hist_usePtGt2020_hypDisamb_usepfMET_usepfJets_useOS_vetoHypMassLt12_requireBTag_sortJetCandidatesbyPt_generalLeptonVeto_createBabyNtuples_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root"){
  setTDRStyle();

  std::cout << "Opening " << FName1.Data() << "\n";
  TFile *f_1         = TFile::Open(FName1.Data());  
  histo1 = AddUpAllMC(f_1,Form("%s_%sj_all",histname.Data(),"all"));  
    
  std::cout << "Opening " << FName2.Data() << "\n";
  TFile *f_2         = TFile::Open(FName2.Data());  
  histo2 = AddUpAllMC(f_2,Form("%s_%sj_all",histname.Data(), "all"));    

 //  std::cout << "Opening " << FName3.Data() << "\n";
//   TFile *f_3         = TFile::Open(FName3.Data());
//   histo3 = AddUpAllMC(f_3,Form("%s_%ij_all",histname.Data(),2));
  
  datahisto = (TH1F*)f_2->Get(Form("data_%s_%sj_all", histname.Data(),"all"));
  std::cout << "Data hist " << datahisto->GetName() << " with entries " << datahisto->GetEntries() << std::endl;
  
  //if several root files loaded in at once, last file loaded doesn't work properly for AddUpAllMC, so do them one at a time.
  
  // double KS1 = histo1->KolmogorovTest(histo2,"UO");
  
  //calculate K-S before rebinning
  double KS1 = histo1->KolmogorovTest(datahisto,"UO");
  double KS2 = histo2->KolmogorovTest(datahisto,"UO");
  // double KS3 = histo3->KolmogorovTest(datahisto,"UO");

  std::cout <<"K-S MC@NLO "<<KS1<< std::endl;  
  std::cout <<"K-S madgraph "<<KS2<< std::endl;  
  //std::cout <<"K-S powheg "<<KS3<< std::endl;
 
  
  datahisto->Rebin(rebin);
  histo1->Rebin(rebin);
  histo2->Rebin(rebin);
  //histo3->Rebin(rebin);


  //double KS1 = histo1->KolmogorovTest(datahisto,"UO");
  //double KS1 = histo1->KolmogorovTest(histo2,"UO");
  //double KS2 = histo2->KolmogorovTest(datahisto,"UO");
  // KS3 = histo3->KolmogorovTest(datahisto,"UO");
  /*
  std::cout <<"K-S MC@NLO "<<KS1<< std::endl;  
  std::cout <<"K-S madgraph "<<KS2<< std::endl;  
  std::cout <<"K-S powheg "<<KS3<< std::endl;
*/

  datahisto->SetMarkerStyle(20);
  datahisto->SetMarkerColor(kBlack);
  datahisto->SetLineColor(kBlack);
  datahisto->GetYaxis()->SetTitle("Events");
  datahisto->GetYaxis()->SetTitleOffset(1.3);	
  datahisto->GetXaxis()->SetTitleOffset(1.0);	

  
  histo1->SetLineColor(kBlue);
  histo1-> SetFillColor(0);
  histo2->SetLineColor(kRed);
  histo2-> SetFillColor(0);
//   histo3->SetLineColor(kGreen+3);
//   histo3-> SetFillColor(0);

  TCanvas *c1 = new TCanvas();
  c1->cd();
  
  datahisto->SetMaximum(1.25*datahisto->GetMaximum());
  if(datahisto->GetMinimum() <0.15 *datahisto->GetMaximum() ) datahisto->SetMinimum(0.);  
  if(datahisto->GetMinimum() > 0.) datahisto->SetMinimum(0.75*datahisto->GetMinimum() );  

  

  if(!drawnorm){
  	datahisto->Draw("Pe");
  	histo1->Draw("histsame");
  	histo2->Draw("histsame");
	//	histo3->Draw("histsame");
  }
  else{
    datahisto->DrawNormalized("Pe");
    histo1->DrawNormalized("histsame");
    histo2->DrawNormalized("histsame");
	//	histo3->DrawNormalized("histsame");
  }

  // TLegend *leg = new TLegend(0.74,0.76,0.92,0.92);
  TLegend *leg = new TLegend(0.4,0.75,0.62,0.88);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  leg->SetTextSize(0.035);
  leg->SetFillStyle(0);
  leg->AddEntry(datahisto, "Data","p");
  // leg->AddEntry(histo3, "powheg","l");
  leg->AddEntry(histo1, "before reweighting","l");
  leg->AddEntry(histo2, "after reweighting","l");

  leg->Draw("same");


  TPaveText *pt1 = new TPaveText(0.20, 0.76, 0.40, 0.91, "brNDC");
  pt1->SetName("pt1name");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  
  TText *blah;
  blah = pt1->AddText("CMS Preliminary, 5.0 fb^{-1} at  #sqrt{s}=7 TeV");
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);  

  // TString KS3_temp = formatFloat(KS3,"%6.2f");
//   KS3_temp.ReplaceAll(" " , "" );
//   KS3_temp = TString("   K-S: ") +  KS3_temp;
//   blah = pt1->AddText(KS3_temp.Data());
//   blah->SetTextSize(0.032);
//   blah->SetTextAlign(11);  
//   blah->SetTextColor(kGreen+3);  

  TString KS1_temp = formatFloat(KS1,"%6.2f");
  KS1_temp.ReplaceAll(" " , "" );
  KS1_temp = TString("   K-S: ") +  KS1_temp;
  blah = pt1->AddText(KS1_temp.Data());
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);  
  blah->SetTextColor(kBlue);  
  /*
  TString KS2_temp = formatFloat(KS2,"%6.2f");
  KS2_temp.ReplaceAll(" " , "" );
  KS2_temp = TString("   K-S: ") +  KS2_temp;
  blah = pt1->AddText(KS2_temp.Data());
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);  
  blah->SetTextColor(kRed);  
  
  */

  
  pt1->Draw();


  c1->Print(Form("mc_comp_%s.pdf", histname.Data()));

  f_1->Close();
  f_2->Close();  
  // f_3->Close();  
  
}

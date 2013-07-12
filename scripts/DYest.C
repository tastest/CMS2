/* 
   Code to do the data driven DY normalization for you.
   Takes in 4 mandatory variables:
   
   dataFName = path with histograms (NJet mainly)
   dataNoMETFName = path with histograms with NO MET cut. This is crucial
   mcFName = path of MC file that has the same cuts applied as dataFName
   scaleMC = lumi to scale MC to
   
   verbose = provides dump of variables used for DY estimate such as
   R_out/in, k_ee, k_mm, etc
*/


#include "TH1F.h"
#include "TString.h"
#include "TCanvas.h"
#include <iostream>
#include "TFile.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "CommonFunctions_DY.C"
#include <fstream>
#include "../tdrStyle.C"
#include <vector>

using namespace std;


TFile *output = new TFile("Routin.root", "RECREATE");
ofstream text; 

TH1F* dyee_hmetOutDYEst_ee;
TH1F* dyee_hmetInDYEst_ee;
TH1F* Ree_vs_met_mc;
TH1F* Ree_vs_met_data;

TH1F* data_hmetInDYEst_ee;
TH1F* data_hmetOutDYEst_ee;

TH1F* dymm_hmetOutDYEst_mm;
TH1F* dymm_hmetInDYEst_mm;
TH1F* Rmm_vs_met_mc;
TH1F* Rmm_vs_met_data;

TH1F* data_hmetInDYEst_mm;
TH1F* data_hmetOutDYEst_mm;

TH1F* data_hmetInDYEst_em;
TH1F* data_hmetOutDYEst_em;

void fill_ratio ( TH1F* & ratio_vs_met, TH1F* hmetOut, TH1F* hmetIn) 
{
  for(int i=1;i<ratio_vs_met->GetNbinsX()+1;i++) {
    float MCo	= GetEntries(hmetOut, i, hmetIn->GetNbinsX()+1);
    float MCoE	= GetTotalError(hmetOut, i, hmetOut->GetNbinsX()+1);
    float MCi	= GetEntries(hmetIn, i, hmetIn->GetNbinsX()+1);
    float MCiE	= GetTotalError(hmetIn, i, hmetIn->GetNbinsX()+1);

    if(MCi == 0) continue;

    float R	= MCo/MCi;
    //float RE	= sqrt(pow(MCoE/MCi,2) + pow(MCo*MCiE/MCi/MCi, 2)/4.); //bug
    float RE	= sqrt(pow(MCoE/MCi,2) + pow(MCo*MCiE/MCi/MCi, 2));
    ratio_vs_met->SetBinContent(i,R);
    ratio_vs_met->SetBinError(i,RE);
  }
/*
  for(int i=1;i<ratio_vs_met->GetNbinsX()+1;i++) {
    float MCo	= GetEntries(hmetOut, i, i);
    float MCoE	= GetTotalError(hmetOut, i, i);
    float MCi	= GetEntries(hmetIn, i, i);
    float MCiE	= GetTotalError(hmetIn, i, i);

    if(MCi == 0) continue;

    float R	= MCo/MCi;
    //float RE	= sqrt(pow(MCoE/MCi,2) + pow(MCo*MCiE/MCi/MCi, 2)/4.); //bug
    float RE	= sqrt(pow(MCoE/MCi,2) + pow(MCo*MCiE/MCi/MCi, 2));
    ratio_vs_met->SetBinContent(i,R);
    ratio_vs_met->SetBinError(i,RE);
  }
*/
}

void fill_ratio_data ( TH1F* & ratio_vs_met, TH1F* hmetOut, TH1F* hmetIn, TH1F* hmetOut_em, TH1F* hmetIn_em, float k, float k_E) 
{
  for(int i=1;i<ratio_vs_met->GetNbinsX()+1;i++) {
    float MCo	= GetEntries(hmetOut, i, hmetIn->GetNbinsX()+1) - 0.5*k*GetEntries(hmetOut_em, i, hmetIn->GetNbinsX()+1) ;
    float MCoE	= sqrt( pow(GetTotalError(hmetOut, i, hmetOut->GetNbinsX()+1),2) + pow(0.5*k*GetTotalError(hmetOut_em, i, hmetOut->GetNbinsX()+1),2) + pow(0.5*k_E*GetEntries(hmetOut_em, i, hmetIn->GetNbinsX()+1),2) );
    float MCi	= GetEntries(hmetIn, i, hmetIn->GetNbinsX()+1) - 0.5*k*GetEntries(hmetIn_em, i, hmetIn->GetNbinsX()+1) ;
    float MCiE	= sqrt( pow(GetTotalError(hmetIn, i, hmetOut->GetNbinsX()+1),2) + pow(0.5*k*GetTotalError(hmetIn_em, i, hmetOut->GetNbinsX()+1),2) + pow(0.5*k_E*GetEntries(hmetIn_em, i, hmetIn->GetNbinsX()+1),2) );    
    
    if(MCi == 0) continue;

    float R	= MCo/MCi;
    //float RE	= sqrt(pow(MCoE/MCi,2) + pow(MCo*MCiE/MCi/MCi, 2)/4.); //bug
    float RE	= sqrt(pow(MCoE/MCi,2) + pow(MCo*MCiE/MCi/MCi, 2)); //ignoring correlation
    ratio_vs_met->SetBinContent(i,R);
    ratio_vs_met->SetBinError(i,RE);
  }
}


double ratio_syst(TH1F* & ratio_vs_met, double ratio, double sfmetcut)
{
  int maxBin = ratio_vs_met->FindBin(sfmetcut);
  
  double dRmax = 0.0;

  for(int i=1;i<maxBin;i++) {
    dRmax = dRmax > TMath::Abs(ratio_vs_met->GetBinContent(i) - ratio) ? dRmax : TMath::Abs(ratio_vs_met->GetBinContent(i) - ratio); 
  }
  return dRmax;
}

void getEstimates( float& n_ee, float& n_mm, 
		   float& n_eeE, float& n_mmE,
		   bool verbose, float scaleMC, 
		   double sfmetcut, int & njet, bool rundata) {

  std::cout << "Get the mininum Bin corresponding to " << sfmetcut << "\n"; 
  double minBin_SF = dyee_hmetInDYEst_ee->FindBin(sfmetcut);
  double maxBin_SF = dyee_hmetInDYEst_ee->GetNbinsX()+1;
  
  std::cout << "minBin_SF " << minBin_SF	 << "\n"; 
  std::cout << "maxBin_SF " << maxBin_SF	 << "\n";
  
  // Do EE
  
  // Calculate the Rout/in for EE in MC
  std::cout << "Calculate the Rout/in for EE in MC...\n";
  float MCo_ee	= GetEntries(dyee_hmetOutDYEst_ee, minBin_SF, maxBin_SF);
  //std::cout << "MCo_ee " << MCo_ee	 << "\n"; 
  float MCo_eeE	= GetTotalError(dyee_hmetOutDYEst_ee, minBin_SF, maxBin_SF);

  float MCi_ee	= GetEntries(dyee_hmetInDYEst_ee, minBin_SF, maxBin_SF);
  float MCi_eeE	= GetTotalError(dyee_hmetInDYEst_ee, minBin_SF, maxBin_SF);

  Ree_vs_met_mc = (TH1F*) dyee_hmetOutDYEst_ee -> Clone(Form("Ree_vs_met_mc_%iJet", njet));
  fill_ratio (Ree_vs_met_mc, dyee_hmetOutDYEst_ee, dyee_hmetInDYEst_ee );

  float MCo_ee_nomet  = GetEntries(dyee_hmetOutDYEst_ee, 0, maxBin_SF);
  float MCo_eeE_nomet = GetTotalError(dyee_hmetOutDYEst_ee, 0, maxBin_SF);

  float MCi_ee_nomet  = GetEntries(dyee_hmetInDYEst_ee, 0, maxBin_SF);
  float MCi_eeE_nomet = GetTotalError(dyee_hmetInDYEst_ee, 0, maxBin_SF);
/*
  float R_ee	= MCo_ee_nomet/MCi_ee_nomet;
  //float R_eeE	= sqrt(pow(MCo_eeE_nomet/MCi_ee_nomet,2) + pow(MCo_ee_nomet*MCi_eeE_nomet/MCi_ee_nomet/MCi_ee_nomet, 2)/4.); //bug
  float R_eeE	= sqrt(pow(MCo_eeE_nomet/MCi_ee_nomet,2) + pow(MCo_ee_nomet*MCi_eeE_nomet/MCi_ee_nomet/MCi_ee_nomet, 2));
  float R_eeE_syst = ratio_syst(Ree_vs_met_mc, R_ee, sfmetcut);
*/  
  //take R_ee after MET cut
  float R_ee	= MCo_ee/MCi_ee;
  float R_eeE	= sqrt(pow(MCo_eeE/MCi_ee,2) + pow(MCo_ee*MCi_eeE/MCi_ee/MCi_ee, 2));
  float R_eeE_syst = ratio_syst(Ree_vs_met_mc, R_ee, sfmetcut);  
  

  // Do the same for data
  std::cout << "Calculate the Rout/in for EE in data...\n";
  float Do_ee	= GetEntries(data_hmetOutDYEst_ee, minBin_SF, maxBin_SF);
  float Do_eeE	= GetTotalError(data_hmetOutDYEst_ee, minBin_SF, maxBin_SF);
  
  //std::cout << "Do_ee " << Do_ee	 << "\n"; 
  
  float Di_ee	= GetEntries(data_hmetInDYEst_ee, minBin_SF, maxBin_SF);
  float Di_eeE	= GetTotalError(data_hmetInDYEst_ee, minBin_SF, maxBin_SF);
  //float Di_eeE	= sqrt(Di_ee);


  // Calculate the Rout/in for MM in MC
  std::cout << "Calculate the Rout/in for MM in MC...\n";
  float MCo_mm	= GetEntries(dymm_hmetOutDYEst_mm, minBin_SF, maxBin_SF);
  float MCo_mmE	= GetTotalError(dymm_hmetOutDYEst_mm, minBin_SF, maxBin_SF);

  float MCi_mm	= GetEntries(dymm_hmetInDYEst_mm, minBin_SF, maxBin_SF);
  float MCi_mmE	= GetTotalError(dymm_hmetInDYEst_mm, minBin_SF, maxBin_SF);
  //float MCi_mmE	= GetTotalError(dymm_hmetOutDYEst_mm, minBin_SF, maxBin_SF); //bug

  Rmm_vs_met_mc = (TH1F*) dymm_hmetOutDYEst_mm -> Clone(Form("Rmm_vs_met_mc_%iJet", njet));
  fill_ratio (Rmm_vs_met_mc, dymm_hmetOutDYEst_mm, dymm_hmetInDYEst_mm );
  
  float MCo_mm_nomet  = GetEntries(dymm_hmetOutDYEst_mm, 0, maxBin_SF);
  float MCo_mmE_nomet = GetTotalError(dymm_hmetOutDYEst_mm, 0, maxBin_SF);

  float MCi_mm_nomet  = GetEntries(dymm_hmetInDYEst_mm, 0, maxBin_SF);
  float MCi_mmE_nomet = GetTotalError(dymm_hmetInDYEst_mm, 0, maxBin_SF);
  /*
  float R_mm	= MCo_mm_nomet/MCi_mm_nomet;
  //float R_mmE	= sqrt(pow(MCo_mmE_nomet/MCi_mm_nomet,2) + pow(MCo_mm_nomet*MCi_mmE_nomet/MCi_mm_nomet/MCi_mm_nomet, 2)/4.); //bug
  float R_mmE	= sqrt(pow(MCo_mmE_nomet/MCi_mm_nomet,2) + pow(MCo_mm_nomet*MCi_mmE_nomet/MCi_mm_nomet/MCi_mm_nomet, 2));
  float R_mmE_syst = ratio_syst(Rmm_vs_met_mc, R_mm, sfmetcut);
  */  
  //take R_mm after MET cut
  float R_mm	= MCo_mm/MCi_mm;
  float R_mmE	= sqrt(pow(MCo_mmE/MCi_mm,2) + pow(MCo_mm*MCi_mmE/MCi_mm/MCi_mm, 2));
  float R_mmE_syst = ratio_syst(Rmm_vs_met_mc, R_mm, sfmetcut);

  // Do the same for data
  std::cout << "Calculate the Rout/in for MM in data...\n";
  float Do_mm	= GetEntries(data_hmetOutDYEst_mm, minBin_SF, maxBin_SF);
  float Do_mmE	= GetTotalError(data_hmetOutDYEst_mm, minBin_SF, maxBin_SF);
  
  //std::cout << "Do_mm " << Do_mm	 << "\n"; 

  float Di_mm	= GetEntries(data_hmetInDYEst_mm, minBin_SF, maxBin_SF);
  float Di_mmE	= GetTotalError(data_hmetInDYEst_mm, minBin_SF, maxBin_SF);
  //float Di_mmE	= sqrt(Di_mm);

  
  // On the OF side, apply the same MET cut 
  // to find the off-peak component to subtract from

  double minBin_OF(minBin_SF), maxBin_OF(maxBin_SF);
  
  std::cout << "minBin_OF " << minBin_OF	 << "\n"; 
  std::cout << "maxBin_OF " << maxBin_OF	 << "\n";
  
  float Di_em		= GetEntries(data_hmetInDYEst_em, minBin_OF, maxBin_OF);
  float Di_emE		= GetTotalError(data_hmetInDYEst_em, minBin_OF, maxBin_OF);
  //float Di_emE	= sqrt(Di_em);
  

  float Do_em		= GetEntries(data_hmetOutDYEst_em, minBin_OF, maxBin_OF);
  float Do_emE		= GetTotalError(data_hmetOutDYEst_em, minBin_OF, maxBin_OF);  
  
  //std::cout << "Do_em " << Do_em	 << "\n"; 
  
    
  // number of SF events without MET cut
  // this is for the data/MC check
  float Di_ee_noMET	= GetEntries(data_hmetInDYEst_ee, 0, maxBin_SF);
  float Di_eeE_noMET	= GetTotalError(data_hmetInDYEst_ee, 0, maxBin_SF);
  //float Di_eeE_noMET	= sqrt(Di_ee_noMET);
  float Di_mm_noMET	= GetEntries(data_hmetInDYEst_mm, 0, maxBin_SF);
  float Di_mmE_noMET	= GetTotalError(data_hmetInDYEst_mm, 0, maxBin_SF);
  //float Di_mmE_noMET	= sqrt(Di_mm_noMET);

/*  
  // we take the e/m reconstruction eff. ratio from the events without MET cut
  float k_ee		= sqrt(Di_ee_noMET/Di_mm_noMET);
  //float k_eeE		= sqrt(pow(Di_eeE_noMET/Di_mm_noMET,2) + pow(Di_mmE_noMET*Di_ee_noMET/Di_mm_noMET/Di_mm_noMET,2)/4.); //bug
  float k_eeE		= k_ee*0.5*sqrt(  pow(Di_eeE_noMET/Di_ee_noMET,2) + pow(Di_mmE_noMET/Di_mm_noMET,2) );
  
  float k_mm		= 1/k_ee;
  //float k_mmE		= sqrt(pow(Di_mmE_noMET/Di_ee_noMET,2) + pow(Di_eeE_noMET*Di_mm_noMET/Di_ee_noMET/Di_ee_noMET,2)/4.); //bug
  float k_mmE		= k_mm*0.5*sqrt(  pow(Di_eeE_noMET/Di_ee_noMET,2) + pow(Di_mmE_noMET/Di_mm_noMET,2) );
*/


  //running on skimmed data, so use ratio with MET cut
  float k_ee		= sqrt(Di_ee/Di_mm);
  float k_eeE		= k_ee*0.5*sqrt(  pow(Di_eeE/Di_ee,2) + pow(Di_mmE/Di_mm,2) );
  
  float k_mm		= 1/k_ee;
  float k_mmE		= k_mm*0.5*sqrt(  pow(Di_eeE/Di_ee,2) + pow(Di_mmE/Di_mm,2) );

  


  Ree_vs_met_data = (TH1F*) data_hmetOutDYEst_ee -> Clone(Form("Ree_vs_met_data_%iJet", njet));
  fill_ratio_data (Ree_vs_met_data, data_hmetOutDYEst_ee, data_hmetInDYEst_ee, data_hmetOutDYEst_em, data_hmetInDYEst_em, k_ee, k_eeE);
  
  Rmm_vs_met_data = (TH1F*) data_hmetOutDYEst_mm -> Clone(Form("Rmm_vs_met_data_%iJet", njet));
  fill_ratio_data (Rmm_vs_met_data, data_hmetOutDYEst_mm, data_hmetInDYEst_mm, data_hmetOutDYEst_em, data_hmetInDYEst_em, k_mm, k_mmE );

  

  float pred_ee	= R_ee*(Di_ee - 0.5*Di_em*k_ee);
  float pred_eeE	= sqrt(pow((Di_ee-0.5*Di_em*k_ee)*R_eeE,2) + pow(R_ee*Di_eeE,2) + pow(0.5*R_ee*k_ee*Di_emE,2) + pow(0.5*R_ee*Di_em*k_eeE,2));
  float pred_eeEc = sqrt(pow((Di_ee-0.5*Di_em*k_ee)*R_eeE,2) + pow(R_ee*Di_eeE*(1-Di_em*k_ee/Di_ee/4),2) + pow(0.5*R_ee*k_ee*Di_emE,2) + pow(k_ee*Di_em*R_ee*Di_mmE/Di_mm/4,2));
  cout<<" pred_eeE excluding and including correlation of k_ee with Di_ee: "<<pred_eeE<<"  "<<pred_eeEc<<endl;
  //float pred_eeE_syst	= sqrt(pow((Di_ee-0.5*Di_em*k_ee)*R_eeE_syst,2) + pow(R_ee*Di_eeE,2) + pow(0.5*R_ee*k_ee*Di_emE,2) + pow(0.5*R_ee*Di_em*k_eeE,2)); //bug
  float pred_eeE_syst	= (Di_ee-0.5*Di_em*k_ee)*R_eeE_syst;

  float pred_mm	= R_mm*(Di_mm - 0.5*Di_em*k_mm);
  float pred_mmE	= sqrt(pow((Di_mm-0.5*Di_em*k_mm)*R_mmE,2) + pow(R_mm*Di_mmE,2) + pow(0.5*R_mm*k_mm*Di_emE,2) + pow(0.5*R_mm*Di_em*k_mmE,2));
  float pred_mmEc = sqrt(pow((Di_mm-0.5*Di_em*k_mm)*R_mmE,2) + pow(R_mm*Di_mmE*(1-Di_em*k_mm/Di_mm/4),2) + pow(0.5*R_mm*k_mm*Di_emE,2) + pow(k_mm*Di_em*R_mm*Di_eeE/Di_ee/4,2));
  cout<<" pred_mmE excluding and including correlation of k_mm with Di_mm: "<<pred_mmE<<"  "<<pred_mmEc<<endl;
  //float pred_mmE_syst	= sqrt(pow((Di_mm-0.5*Di_em*k_mm)*R_mmE_syst,2) + pow(R_mm*Di_mmE,2) + pow(0.5*R_mm*k_mm*Di_emE,2) + pow(0.5*R_mm*Di_em*k_mmE,2)); //bug
  float pred_mmE_syst	= (Di_mm-0.5*Di_em*k_mm)*R_mmE_syst;

  if(verbose) {
    cout << "Do_ee: " << Do_ee << "+/-" << Do_eeE << endl;
    cout << "Do_mm: " << Do_mm << "+/-" << Do_mmE << endl;
    cout << "Do_em: " << Do_em << "+/-" << Do_emE << endl;
    cout << "D_em: " << Do_em + Di_em << "+/-" << sqrt(Do_emE*Do_emE + Di_emE*Di_emE) << endl;
    cout << endl;
    
    cout << "Di_ee: " << Di_ee << "+/-" << Di_eeE << endl;
    cout << "R_ee: " << R_ee << "+/-" << R_eeE << endl;
    cout << "k_ee: " << k_ee << "+/-" << k_eeE << endl;
    cout << "Di_em: " << Di_em << "+/-" << Di_emE << endl;
    cout<<"pred_ee: "<<pred_ee<<"+/-"<<pred_eeE<<"+/-"<<pred_eeE_syst<<endl;
        
    cout << "Di_mm: " << Di_mm << "+/-" << Di_mmE << endl;
    cout << "R_mm: " << R_mm << "+/-" << R_mmE << endl;
    cout << "k_mm: " << k_mm << "+/-" << k_mmE << endl;
    cout << "Di_em: " << Di_em << "+/-" << Di_emE << endl;
    cout<<"pred_mm: "<<pred_mm<<"+/-"<<pred_mmE<<"+/-"<<pred_mmE_syst<<endl;
    
    double scale_ee = GetEntries(dyee_hmetOutDYEst_ee, 0, dyee_hmetOutDYEst_ee->GetNbinsX()+1)/dyee_hmetOutDYEst_ee->GetEntries(); // get the weight of the histogram
    double scale_mm = GetEntries(dymm_hmetOutDYEst_mm, 0,  dymm_hmetOutDYEst_mm->GetNbinsX()+1)/dymm_hmetOutDYEst_mm->GetEntries(); // get the weight of the histogram

    cout << "MCi_ee: " << MCi_ee << "+/-" << MCi_eeE << "; event count = "<< Form("%4.1f", MCi_ee/scale_ee)<< endl;
    cout << "MCo_ee: " << MCo_ee << "+/-" << MCo_eeE << "; event count = "<< Form("%4.1f", MCo_ee/scale_ee)<< endl;

    cout << "MCi_mm: " << MCi_mm << "+/-" << MCi_mmE << "; event count = "<< Form("%4.1f", MCi_mm/scale_mm)<< endl;
    cout << "MCo_mm: " << MCo_mm << "+/-" << MCo_mmE << "; event count = "<< Form("%4.1f", MCo_mm/scale_mm)<< endl;
  }
  cout << endl;
  text << "\\begin{tabular}{l |  c    c}" << endl; 
  text << "\\hline" << endl;
  text << "Sample & ee & $\\mu\\mu$  \\\\ \\hline" << endl;  
  text << "$N_{in}$ in data & " << Form("%4.2f %s %4.2f", Di_ee, " $\\pm$", Di_eeE) << " & "
       << Form("%4.2f %s %4.2f", Di_mm, " $\\pm$", Di_mmE) << " \\\\ " << endl;

  text << "$N_{in}$ in $e\\mu$ data & \\multicolumn{2}{c}{ " << Form("%4.2f %s %4.2f", Di_em, " $\\pm$", Di_emE) 
       << "  } \\\\ " << endl;

  text << "$k$ & " << Form("%4.2f %s %4.2f", k_ee, " $\\pm$", k_eeE) << " & "
       << Form("%4.2f %s %4.2f", k_mm, " $\\pm$", k_mmE) << " \\\\ " << endl;
       
  text << "$R_{out/in}$ & " << Form("%4.2f", R_ee) 
       << " $\\pm$ " << Form("%4.2f", R_eeE); 
  if (sfmetcut != 0) 
    text << " $\\pm$ " << Form("%4.2f", R_eeE_syst) ;
  text  << " & "  << Form("%4.2f", R_mm) 
       << " $\\pm$ " << Form("%4.2f", R_mmE) ;
  if(sfmetcut != 0) 
    text << " $\\pm$ " << Form("%4.2f", R_mmE_syst);
  //text << "   \\\\" <<endl;


  text   << "\\\\ \\hline" << endl;   

  if(MCi_ee > 0) {
    text << "DY Estimate (preselection) & " << Form("%4.2f", pred_ee) 
	 << " $\\pm$ " << Form("%4.2f", pred_eeE);
    if(sfmetcut != 0) 
      text << " $\\pm$  " << Form("%4.2f", pred_eeE_syst);
    text << " & ";
  }
  else 
    text << "DY Estimate (preselection) & " << 0 << " $\\pm$ " << 0 << " & ";
  if(MCi_mm > 0) {
    text << Form("%4.2f", pred_mm) << " $\\pm$ " << Form("%4.2f", pred_mmE);
    if (sfmetcut != 0)
      text  << " $\\pm$ " << Form("%4.2f", pred_mmE_syst);
  }
  else 
    text << 0 << " $\\pm$ " << 0;
  text   << " \\\\ " << endl;
    
  if(!rundata){
  text << "DY MC event count & " << Form("%4.2f", scaleMC*GetEntries(dyee_hmetOutDYEst_ee, minBin_SF, maxBin_SF)) 
       << " $\\pm$ " << Form("%4.2f", scaleMC*GetTotalError(dyee_hmetOutDYEst_ee, minBin_SF, maxBin_SF)) << " & ";
  text << Form("%4.2f", scaleMC*GetEntries(dymm_hmetOutDYEst_mm, minBin_SF, maxBin_SF)) 
       << " $\\pm$ " << Form("%4.2f", scaleMC*GetTotalError(dymm_hmetOutDYEst_mm, minBin_SF, maxBin_SF)) << "\\\\ "<<endl;
  }
  
  text   << " \\hline" << endl; 
 
  /*
  //text << "N(out) in data & " << Form("%4.2f", Do_ee) << " & " << Form("%4.2f", Do_mm) 
  //     << "\\\\ \\hline" << endl;
       
  text << "N(out) in data & " << Form("%4.2f %s %4.2f", Do_ee, " $\\pm$", Do_eeE) << " & "
       << Form("%4.2f %s %4.2f", Do_mm, " $\\pm$", Do_mmE) << " \\\\ " << endl;       
  text << "N(out) in $e\\mu$ data & " << Form("%4.2f %s %4.2f", Do_em, " $\\pm$", Do_emE) 
       << " \\\\ \\hline" << endl;       
  */

/*
  //p_T>30   
  double S_ee = 0.072644913;
  double S_eeE = 0.011646098;  
  double S_mm = 0.062299242;
  double S_mmE = 0.009000063;

  //p_T>50
  //double S_ee = 0.135533477;
  //double S_eeE = 0.027417464;  
  //double S_mm = 0.099232846;
  //double S_mmE = 0.019675574;
*/


  //MET>50   
  double S_ee = 0.108423523;
  double S_eeE = 0.039208951;  
  double S_mm = 0.14623613;
  double S_mmE = 0.035346489;




text << "$m_{lb}^{min}$ mass cut efficiency for DY MC & " << Form("%4.2f %s %4.2f", S_ee, " $\\pm$", S_eeE) << " & "
       << Form("%4.2f %s %4.2f", S_mm, " $\\pm$", S_mmE) << " \\\\ " << endl;   
  
  if(MCi_ee > 0) {
    text << "DY Estimate (full selection) & " << Form("%4.2f", pred_ee*S_ee) 
	 << " $\\pm$ " << Form("%4.2f", sqrt( pow( (pred_eeE/pred_ee) ,2) + pow( (S_eeE/S_ee) ,2) )*fabs(pred_ee)*S_ee  );
    if(sfmetcut != 0) 
      text << " $\\pm$  " << Form("%4.2f", pred_eeE_syst*S_ee);
    text << " & ";
  }
  else 
    text << "DY Estimate (full selection) & " << 0 << " $\\pm$ " << 0 << " & ";
  if(MCi_mm > 0) {
    text << Form("%4.2f", pred_mm*S_mm) << " $\\pm$ " << Form("%4.2f", sqrt( pow( (pred_mmE/pred_mm) ,2) + pow( (S_mmE/S_mm) ,2) )*fabs(pred_mm)*S_mm );
    if (sfmetcut != 0)
      text  << " $\\pm$ " << Form("%4.2f", pred_mmE_syst*S_mm);
  }
  else 
    text << 0 << " $\\pm$ " << 0;

  text   << "\\\\ \\hline" << endl;
       
  text << "\\end{tabular}" << endl;
  
  n_ee = pred_ee;
  n_eeE = pred_eeE;
  n_mm = pred_mm;
  n_mmE = pred_mmE;
  
  //cout<<pred_ee<<endl;
  //cout<<pred_mm<<endl;
 
  cout<<"upper limit ee preselection "<< pred_ee + sqrt(pred_eeE*pred_eeE + pred_eeE_syst*pred_eeE_syst)<<endl;
  cout<<"upper limit mm preselection "<< pred_mm + sqrt(pred_mmE*pred_mmE + pred_mmE_syst*pred_mmE_syst)<<endl;
 
  cout<<"DY Estimate ee (full selection): "<<pred_ee*S_ee<<endl;
  cout<<"DY Estimate mm (full selection): "<<pred_mm*S_mm<<endl;
  
  cout<<"upper limit ee full selection "<< pred_ee*S_ee + sqrt(  ( pow( (pred_eeE/pred_ee) ,2) + pow( (S_eeE/S_ee) ,2) )*pred_ee*pred_ee*S_ee*S_ee + pred_eeE_syst*pred_eeE_syst*S_ee*S_ee )<<endl;
  cout<<"upper limit mm full selection "<< pred_mm*S_mm + sqrt(  ( pow( (pred_mmE/pred_mm) ,2) + pow( (S_mmE/S_mm) ,2) )*pred_mm*pred_mm*S_mm*S_mm + pred_mmE_syst*pred_mmE_syst*S_mm*S_mm )<<endl;
  cout << endl << endl;
  
  float eeplusmumu_pre = R_ee*Di_ee + R_mm*Di_mm - 0.5*Di_em*(R_ee*k_ee + R_mm*k_mm );
  //can't just sum errors in quadrature due to correlations (both predictions depend on Di_em and k)
  float eeplusmumu_preE = sqrt( pow(R_eeE*(Di_ee-0.5*Di_em*k_ee),2) + pow(Di_eeE*R_ee,2)  +  pow(R_mmE*(Di_mm-0.5*Di_em*k_mm),2) + pow(Di_mmE*R_mm,2)
                          + pow(0.5*Di_emE*(R_ee*k_ee + R_mm*k_mm ),2)  + pow(0.5*Di_em*(   (R_mm/k_ee/k_ee-R_ee)*k_eeE   ),2) ) ;
                          
  //this is correct when k_ee=sqrt(Di_ee/Di_mm), but eeplusmumu_preE gives practically the same result
  //float eeplusmumu_preE_new =   sqrt(  pow(R_eeE*(Di_ee-0.5*Di_em*k_ee),2) + pow(Di_eeE*(R_ee*(1-Di_em/4/k_ee/Di_mm)+Di_em*Di_mm*R_mm*k_ee/4/Di_ee/Di_ee),2) + pow(Di_emE*0.5*(k_ee*R_ee + k_mm*R_mm),2) + pow(Di_mmE*(R_mm*(1-Di_em/4/Di_ee/k_mm)+Di_ee*Di_em*R_ee*k_mm/4/Di_mm/Di_mm),2) + pow(R_mmE*(Di_mm-0.5*Di_em*k_mm),2)  );
                          
                                                    

  float eeplusmumu_full = R_ee*Di_ee*S_ee + R_mm*Di_mm*S_mm - 0.5*Di_em*(R_ee*k_ee*S_ee + R_mm*k_mm*S_mm ) ;
  //can't just sum errors in quadrature due to correlations
  float eeplusmumu_fullE = sqrt( pow(R_eeE*(Di_ee-0.5*Di_em*k_ee)*S_ee,2) + pow(Di_eeE*R_ee*S_ee,2)  +  pow(R_mmE*(Di_mm-0.5*Di_em*k_mm)*S_mm,2) + pow(Di_mmE*R_mm*S_mm,2)
                          + pow(0.5*Di_emE*(R_ee*k_ee*S_ee + R_mm*k_mm*S_mm ),2)  + pow(0.5*Di_em*(   (R_mm*S_mm/k_ee/k_ee-R_ee*S_ee)*k_eeE   ),2)  
                          + pow(S_eeE*(R_ee*Di_ee -  0.5*Di_em*R_ee*k_ee),2) + pow(S_mmE*(R_mm*Di_mm - 0.5*Di_em*R_mm*k_mm),2)  );

  
    
  //cout<<"sum of ee and mumu (preselection) "<<eeplusmumu_pre<<"   "<<pred_ee + pred_mm<<endl;
  //cout<<"stat uncertainty: "<<eeplusmumu_preE<<" stat uncertainty new: "<<eeplusmumu_preE_new<<"   simple sum in quadrature: "<<sqrt( pred_eeE*pred_eeE + pred_mmE*pred_mmE)<<endl;
  //
  //cout<<"sum of ee and mumu (full selection) "<<eeplusmumu_full<<"   "<< pred_ee*S_ee + pred_mm*S_mm<<endl;
  //cout<<"stat uncertainty: "<<eeplusmumu_fullE<<"   simple sum in quadrature: "<<sqrt( pow(pred_eeE*S_ee,2) + pow(pred_ee*S_eeE,2) + pow(pred_mmE*S_mm,2) + pow(pred_mm*S_mmE,2) )<<endl;

    
  cout<<"sum of ee and mumu (preselection) "<<pred_ee + pred_mm<<endl;
  cout<<"stat uncertainty: "<<eeplusmumu_preE<<"   syst uncertainty: "<<sqrt( pred_eeE_syst*pred_eeE_syst + pred_mmE_syst*pred_mmE_syst)<<"    total uncertainty: "<<sqrt( pow(eeplusmumu_preE,2) + pow(pred_eeE_syst,2) + pow(pred_mmE_syst,2) ) <<endl;
  
  cout<<"sum of ee and mumu (full selection) "<<pred_ee*S_ee + pred_mm*S_mm<<endl;
  cout<<"stat uncertainty: "<<eeplusmumu_fullE<<"   syst uncertainty: "<<sqrt( pow(pred_eeE_syst*S_ee,2) + pow(pred_mmE_syst*S_mm,2) )<<"    total uncertainty: "<<sqrt( pow(eeplusmumu_fullE,2) + pow(pred_eeE_syst*S_ee,2) + pow(pred_mmE_syst*S_mm,2) ) <<endl<<endl;
  
  
  setStyle(Ree_vs_met_mc, 1, kBlue, 20, "MET cut (GeV)", "R_{out/in}", 30, 50);
  setStyle(Rmm_vs_met_mc, 1, kRed, 20, "MET cut (GeV)", "R_{out/in}", 30, 50);

  setStyle(Ree_vs_met_data, 1, kBlue, 20, "MET cut (GeV)", "R_{out/in}", 30, 50);
  setStyle(Rmm_vs_met_data, 1, kRed, 20, "MET cut (GeV)", "R_{out/in}", 30, 50);

  Ree_vs_met_mc->GetYaxis()->SetRangeUser(0,0.5);
  Rmm_vs_met_mc->GetYaxis()->SetRangeUser(0,0.5);

  Ree_vs_met_data->GetYaxis()->SetRangeUser(-1,2.0);
  Rmm_vs_met_data->GetYaxis()->SetRangeUser(-1,2.0);


  TCanvas *c1 = new TCanvas();
  c1->cd();  

  TLegend *leg = new TLegend(0.2, 0.75, 0.4, 0.90, "", "brNDC");
  leg->SetFillColor(0);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetShadowColor(0);
  leg->AddEntry(Ree_vs_met_mc, "ee", "lp");
  leg->AddEntry(Rmm_vs_met_mc, "#mu#mu", "lp");
  Ree_vs_met_mc->Draw("h");
  Rmm_vs_met_mc->Draw("same e1");
  leg->Draw();
  c1->Print(Form("Routin_mc_%iJet.eps", njet));
  c1->Print(Form("Routin_mc_%iJet.png", njet));
  c1->Print(Form("Routin_mc_%iJet.pdf", njet));



  TCanvas *c2 = new TCanvas();
  c2->cd();  
  TLegend *leg2 = new TLegend(0.2, 0.75, 0.4, 0.90, "", "brNDC");
  leg2->SetFillColor(0);
  leg2->AddEntry(Ree_vs_met_data, "ee", "lp");
  leg2->AddEntry(Rmm_vs_met_data, "#mu#mu", "lp");
  Ree_vs_met_data->Draw("h");
  Rmm_vs_met_data->Draw("same e1");
  leg2->Draw();
  c2->Print(Form("Routin_data_%iJet.eps", njet));
  c2->Print(Form("Routin_data_%iJet.png", njet));
  c2->Print(Form("Routin_data_%iJet.pdf", njet));
  




  TCanvas *c3 = new TCanvas();
  c3->Divide(4,3);
  //gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0010);
  gStyle->SetStatW(0.6);    
  gStyle->SetStatH(0.002);   
  gStyle->SetStatY(1.05);
  c3->cd(1);
  dyee_hmetOutDYEst_ee->Draw();
  c3->cd(2);
  dyee_hmetInDYEst_ee->Draw();
  c3->cd(3);
  dymm_hmetOutDYEst_mm->Draw();
  c3->cd(4);
  dymm_hmetInDYEst_mm->Draw();
  c3->cd(5);
  data_hmetOutDYEst_ee->Draw();
  c3->cd(6);
  data_hmetInDYEst_ee->Draw();
  c3->cd(7);
  data_hmetOutDYEst_mm->Draw();
  c3->cd(8);
  data_hmetInDYEst_mm->Draw();
  c3->cd(10);
  data_hmetOutDYEst_em->Draw();
  c3->cd(11);
  data_hmetInDYEst_em->Draw();

  
  
  c3->Print(Form("Histograms_%iJet.eps", njet));
  c3->Print(Form("Histograms_%iJet.png", njet));
  c3->Print(Form("Histograms_%iJet.pdf", njet));
  



  output->cd();
  Ree_vs_met_mc->Write();
  Rmm_vs_met_mc->Write();

  Ree_vs_met_data->Write();
  Rmm_vs_met_data->Write();
  // output->Close();
  
}

TH1F* AddUpAllMC(TFile* f_mc, char* histname) {
	std::vector<TString> mcsamplenames;
	mcsamplenames.push_back("ttdil");
	mcsamplenames.push_back("ttotr");
	mcsamplenames.push_back("wjets");
	mcsamplenames.push_back("DYee");
	mcsamplenames.push_back("DYmm");
	mcsamplenames.push_back("DYtautau");
	mcsamplenames.push_back("VV");
	mcsamplenames.push_back("tw");
	
	// std::cout << Form("%s_%s",mcsamplenames[0].Data(),histname) << std::endl;

	TH1F *clone = (TH1F*)f_mc->Get(Form("%s_%s",mcsamplenames[0].Data(),histname))->Clone();
	// std::cout << "adding hist " << clone->GetName() << " with entries " << clone->GetEntries() << std::endl;	
	for (unsigned int i = 1; i < mcsamplenames.size(); ++i) {
		TH1F *hist = (TH1F*)f_mc->Get(Form("%s_%s",mcsamplenames[i].Data(),histname));
		// std::cout << "adding hist " << hist->GetName() << " with entries " << hist->GetEntries() << std::endl;
		clone->Add(hist);
	}
	
	return clone;	
}



void DYest(TString dataFName, TString mcFName, bool rundata = true, float METcut = 30., float scaleMC=1., bool verbose = true){
  setTDRStyle(); 

  std::cout << "Opening " << dataFName.Data() << "\n";
  TFile *f_data		= TFile::Open(dataFName.Data());
  std::cout << "Opening " << mcFName.Data() << "\n";
  TFile *f_mc		= TFile::Open(mcFName.Data());


  
  text.open("DYest.tex");
  text << "\\documentclass{article}" << endl;
  text << "\\usepackage{times}" << endl;
  text << "\\usepackage{epsfig}" << endl;
  text << "\\begin{document}" << endl;


  for(int njet=2;njet<3;njet++) {
  	
  	
  	if(rundata) {  	
    std::cout << "Opening " << Form("dyee_hmetOutDYEst_%ij_ee",njet) << "\n";

    dyee_hmetOutDYEst_ee	= (TH1F*)f_mc->Get(Form("DYee_hmetOutDYEst_%ij_ee",njet));
    std::cout << "dyee_hmetOutDYEst_ee -> GetEntries() = " << dyee_hmetOutDYEst_ee -> GetEntries() << endl;
    dyee_hmetInDYEst_ee	        = (TH1F*)f_mc->Get(Form("DYee_hmetInDYEst_%ij_ee",njet));
    data_hmetOutDYEst_ee	= (TH1F*)f_data->Get(Form("data_hmetOutDYEst_%ij_ee",njet));
    data_hmetInDYEst_ee  	= (TH1F*)f_data->Get(Form("data_hmetInDYEst_%ij_ee",njet));

    dymm_hmetOutDYEst_mm	= (TH1F*)f_mc->Get(Form("DYmm_hmetOutDYEst_%ij_mm",njet));
    dymm_hmetInDYEst_mm	        = (TH1F*)f_mc->Get(Form("DYmm_hmetInDYEst_%ij_mm",njet));
    data_hmetOutDYEst_mm	= (TH1F*)f_data->Get(Form("data_hmetOutDYEst_%ij_mm",njet));
    data_hmetInDYEst_mm	        = (TH1F*)f_data->Get(Form("data_hmetInDYEst_%ij_mm",njet));
    
    //get the emu data hists
    data_hmetInDYEst_em	        = (TH1F*)f_data->Get(Form("data_hmetInDYEst_%ij_em",njet));
    data_hmetOutDYEst_em	        = (TH1F*)f_data->Get(Form("data_hmetOutDYEst_%ij_em",njet));
    }
    
   	else {
    std::cout << "Opening " << Form("dyee_hmetOutDYEst_%ij_ee",njet) << "\n";

    dyee_hmetOutDYEst_ee	= (TH1F*)f_mc->Get(Form("DYee_hmetOutDYEst_%ij_ee",njet));
    std::cout << "dyee_hmetOutDYEst_ee -> GetEntries() = " << dyee_hmetOutDYEst_ee -> GetEntries() << endl;
    dyee_hmetInDYEst_ee	        = (TH1F*)f_mc->Get(Form("DYee_hmetInDYEst_%ij_ee",njet));
    data_hmetOutDYEst_ee	= AddUpAllMC(f_mc,Form("hmetOutDYEst_%ij_ee",njet));
    data_hmetInDYEst_ee  	= AddUpAllMC(f_mc,Form("hmetInDYEst_%ij_ee",njet));

    dymm_hmetOutDYEst_mm	= (TH1F*)f_mc->Get(Form("DYmm_hmetOutDYEst_%ij_mm",njet));
    dymm_hmetInDYEst_mm	        = (TH1F*)f_mc->Get(Form("DYmm_hmetInDYEst_%ij_mm",njet));
    data_hmetOutDYEst_mm	= AddUpAllMC(f_mc,Form("hmetOutDYEst_%ij_mm",njet));
    data_hmetInDYEst_mm	        = AddUpAllMC(f_mc,Form("hmetInDYEst_%ij_mm",njet));
    
    //get the emu data hists
    data_hmetInDYEst_em	        = AddUpAllMC(f_mc,Form("hmetInDYEst_%ij_em",njet));
    data_hmetOutDYEst_em	        = AddUpAllMC(f_mc,Form("hmetOutDYEst_%ij_em",njet));
    }



    float n_ee = 0;
    float n_eeE =0; 
    float n_mm = 0;
    float n_mmE =0;
           
    if(verbose) {
      cout << endl << "//-----------------------------------------------------------------------//" << endl << endl;
      cout << "Doing the " << njet << " jet Bin:" << endl;
    }
                
    // Writing the table without MET cut
    std::cout << "Writing the table without MET cut...\n";
    text << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
    text << "\\begin{table}" << endl;	
  	if(rundata) {  	
    text << "\\caption{Drell-Yan estimation with 30 GeV MET Cut in the " << njet << " Jet bin. }"<<endl;
  	}
  	else {
  	text << "\\caption{Drell-Yan estimation MC closure test with 30 GeV MET Cut in the " << njet << " Jet bin. }"<<endl;  	}
    getEstimates(n_ee, n_mm, n_eeE, n_mmE, verbose, scaleMC, 0, njet, rundata);
    text << "\\end{table}" << endl;
    
    std::cout << "Writing the table with 50 GeV MET cut...\n";
    text << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
    text << "\\begin{table}" << endl;  	
  	if(rundata) {  	
    text << "\\caption{Drell-Yan estimation with 50 GeV MET Cut in the " << njet << " Jet bin. }"<<endl;
  	}
  	else {
  	text << "\\caption{Drell-Yan estimation closure test with 50 GeV MET Cut in the " << njet << " Jet bin. }"<<endl;  	}
    getEstimates(n_ee, n_mm, n_eeE, n_mmE, verbose, scaleMC, METcut, njet, rundata);
    text << "\\end{table}" << endl;    

    text << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
    text << "\\begin{figure}[htb]"<<endl;
    text << "\\begin{center}"<<endl;
    text << "\\begin{tabular}{c}"<<endl;
    text << "\\epsfig{figure=Routin_mc_" << njet<< "Jet.pdf, width=2.5in}"<<endl; 
    text << "\\epsfig{figure=Routin_data_"<< njet << "Jet.pdf, width=2.5in}"<<endl; 
    text << "\\end{tabular}"<<endl;
    text << "\\caption{$R_{out/in}$ dependence on the MET cut in the " << njet << " Jet bin. }"<<endl;
    text << "\\end{center}"<<endl;
    text << "\\end{figure}"<<endl;
    text << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
    
  }
  
  
  text << "\\end{document}" << endl;

  f_data->Close();
  f_mc->Close();
  
}






#include "plotResults.h"
#include "../../Tools/HistogramUtilities.h"
#include "../../Tools/AllDataSources.h"
#include "../../Tools/Utilities.h"
#include "../../Tools/histtools.cc"

#include "TROOT.h"

#include "TGraphAsymmErrors.h"
#include "TArrow.h"
#include "TLatex.h"

// for 2_2_X
const static sources_t theMCSignal =  (1<<H_WJETS);

const static sources_t theMCBackground =  (1<<H_QCD30);

const static sources_t theMCSources = theMCSignal | theMCBackground;

const static sources_t theSources = theMCSignal | theMCBackground;

const static sources_t theDataSources = (1<<H_WHUNT);



TArrow *getArrow(THStack *st, TString det, float cutValEB, float cutValEE, float max)
{

    float cutVal; 
    if (det == "eb") cutVal = cutValEB;
    else cutVal = cutValEE;
    if (max == -1.0) max = st->GetMaximum()/1.2;
    TArrow *arr_cut = new TArrow(cutVal, max, cutVal, 0, 0.08, "|>");
    arr_cut->SetAngle(25.0);
    arr_cut->SetLineWidth(3.0);
    arr_cut->SetLineColor(kRed);
    arr_cut->SetFillColor(kRed);
    arr_cut->SetLineStyle(kDashed);
    arr_cut->SetFillStyle(3001);
    return arr_cut;
}


void plotStack(HistogramUtilities &h1, TString name, TString titleX, TString saveName, TString det, int rebin, float cutValEB, float cutValEE)
{

    THStack *st = h1.getStack(theSources, name, "", det, rebin);
    TLegend *lg_all = h1.getLegend(theSources, name, "", det);
    lg_all->SetX1(0.55);
    lg_all->SetX2(0.70);
    lg_all->SetY1(0.6);
    lg_all->SetY2(0.8);

    lg_all->SetTextSize(0.04);

    TCanvas *c1 = new TCanvas();
    c1->cd();
    st->Draw();
    st->GetXaxis()->SetTitle(titleX); 
    lg_all->Draw();
    if (det == "ee" && cutValEE != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    if (det == "eb" && cutValEB != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    if (det == "all" && cutValEE != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    Utilities::saveCanvas(c1, "results/" + saveName  + "_lin_" + name + "_" + det);

    c1->SetLogy();
    st->SetMinimum(0.10);
    st->Draw();	
    lg_all->Draw();
    if (det == "ee" && cutValEE != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    if (det == "eb" && cutValEB != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    if (det == "all" && cutValEE != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    Utilities::saveCanvas(c1, "results/" + saveName  + "_log_" + name + "_" + det);

    delete c1;
    delete st;
    delete lg_all;

}


void plotDataRefOverlayStack(HistogramUtilities &hData, HistogramUtilities &hRef, 
        TString name, TString titleX, TString saveName, TString det, int rebin,  float cutValEB, float cutValEE)
{

    TLatex *   tex = new TLatex(0.55,0.83,"MC Norm (1 nb^{-1})");
    tex->SetNDC();
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    TLatex *   tex2 = new TLatex(0.55,0.88,"CMS Preliminary 2010");
    tex2->SetNDC();
    tex2->SetTextSize(0.04);
    tex2->SetLineWidth(2);


    THStack *stRef = hRef.getStack(theMCSources, name, "", det, rebin);
    TLegend *lg_all = hRef.getLegend(theMCSources, name, "", det);
    lg_all->SetX1(0.55);
    lg_all->SetX2(0.70);
    lg_all->SetY1(0.6);
    lg_all->SetY2(0.8);
    lg_all->SetTextSize(0.04);

    TH1F *h1Data =  hData.getHistogram(theDataSources, name, "", det, rebin);
    h1Data->GetXaxis()->SetTitle("");
    h1Data->GetYaxis()->SetTitle("");
    h1Data->GetXaxis()->SetLabelSize(0);
    h1Data->GetYaxis()->SetLabelSize(0);
    h1Data->SetMarkerStyle(20);
    h1Data->SetMarkerSize(1.5);
    h1Data->SetLineWidth(2.0);
    h1Data->SetMarkerColor(kRed);
    h1Data->SetLineColor(kRed);
    lg_all->AddEntry(h1Data,"Data", "lep"); 

    // ideal y range should be at least 2 sigma greater than highest data point
    float idealMaximumData = h1Data->GetMaximum() + sqrt(h1Data->GetMaximum())*2;
    if(stRef->GetMaximum() < idealMaximumData ) stRef->SetMaximum(idealMaximumData);

    TCanvas *c1 = new TCanvas();
    c1->cd();
    stRef->Draw();
    h1Data->Draw("samee1");
    stRef->GetXaxis()->SetTitle(titleX); 
    lg_all->Draw();
    tex->Draw();
    tex2->Draw();
    if (det == "ee" && cutValEE != -1) getArrow(stRef, det, cutValEB, cutValEE, idealMaximumData/1.2)->Draw();
    if (det == "eb" && cutValEB != -1) getArrow(stRef, det, cutValEB, cutValEE, idealMaximumData/1.2)->Draw();
    if (det == "all" && cutValEE != -1) getArrow(stRef, det, cutValEB, cutValEE, idealMaximumData/1.2)->Draw();
    Utilities::saveCanvas(c1, "results/" + saveName  + "_lin_" + name + "_" + det);

    c1->SetLogy();
    stRef->SetMinimum(0.10);
    stRef->Draw();	
    h1Data->Draw("samee1");
    lg_all->Draw();
    tex->Draw();
    tex2->Draw();
    if (det == "ee" && cutValEE != -1) getArrow(stRef, det, cutValEB, cutValEE, idealMaximumData/1.2)->Draw();
    if (det == "eb" && cutValEB != -1) getArrow(stRef, det, cutValEB, cutValEE, idealMaximumData/1.2)->Draw();
    if (det == "all" && cutValEE != -1) getArrow(stRef, det, cutValEB, cutValEE, idealMaximumData/1.2)->Draw();
    Utilities::saveCanvas(c1, "results/" + saveName  + "_log_" + name + "_" + det);

    delete tex;
    delete tex2;
    delete c1;
    delete stRef;
    delete lg_all;

}

void plotElectronIDStack(HistogramUtilities &hData, HistogramUtilities &hRef, TString name, TString titleX, TString saveName, TString det, int rebin)
{

    TString selected_name = "ele_selected_" + name;
    TString antiselected_name = "ele_antiselected_" + name;

    TLatex *   texSelected = new TLatex(0.55,0.88,"MET > 20.0 GeV");
    texSelected->SetNDC();
    texSelected->SetTextSize(0.04);
    texSelected->SetLineWidth(2);

    TLatex *   texAntiselected = new TLatex(0.55,0.88,"MET < 15.0 GeV");
    texAntiselected->SetNDC();
    texAntiselected->SetTextSize(0.04);
    texAntiselected->SetLineWidth(2);

    TLatex *   texLumi = new TLatex(0.55,0.83,"MC Norm (1 nb^{-1})");
    texLumi->SetNDC();
    texLumi->SetTextSize(0.04);
    texLumi->SetLineWidth(2);

    // get selected and antiselected reference mc
    THStack *stRefSelected = hRef.getStack(theMCSources, selected_name, "", det, rebin);
    THStack *stRefAntiselected = hRef.getStack(theMCSources, antiselected_name, "", det, rebin);

    // get selected data
    TH1F *h1DataSelected =  hData.getHistogram(theDataSources, selected_name, "", det, rebin);
    h1DataSelected->GetXaxis()->SetTitle("");
    h1DataSelected->GetYaxis()->SetTitle("");
    h1DataSelected->GetXaxis()->SetLabelSize(0);
    h1DataSelected->GetYaxis()->SetLabelSize(0);
    h1DataSelected->SetMarkerStyle(20);
    h1DataSelected->SetMarkerSize(1.5);
    h1DataSelected->SetLineWidth(2.0);
    h1DataSelected->SetMarkerColor(kRed);
    h1DataSelected->SetLineColor(kRed);

    // get antiselected data
    TH1F *h1DataAntiselected =  hData.getHistogram(theDataSources, antiselected_name, "", det, rebin);
    h1DataAntiselected->GetXaxis()->SetTitle("");
    h1DataAntiselected->GetYaxis()->SetTitle("");
    h1DataAntiselected->GetXaxis()->SetLabelSize(0);
    h1DataAntiselected->GetYaxis()->SetLabelSize(0);
    h1DataAntiselected->SetMarkerStyle(20);
    h1DataAntiselected->SetMarkerSize(1.5);
    h1DataAntiselected->SetLineWidth(2.0);
    h1DataAntiselected->SetMarkerColor(kRed);
    h1DataAntiselected->SetLineColor(kRed);

    TLegend *lg_all = hRef.getLegend(theMCSources, antiselected_name, "", det);
    lg_all->SetX1(0.55);
    lg_all->SetX2(0.70);
    lg_all->SetY1(0.6);
    lg_all->SetY2(0.8);
    lg_all->SetTextSize(0.04);
    lg_all->AddEntry(h1DataSelected, "Data", "lep");

    TCanvas *c1 = new TCanvas();
    c1->cd();

    // draw selected
    // ideal y range should be at least 2 sigma greater than highest data point
    float idealMaximumData = h1DataSelected->GetMaximum() + sqrt(h1DataSelected->GetMaximum())*2;
    if(stRefSelected->GetMaximum() < idealMaximumData ) stRefSelected->SetMaximum(idealMaximumData);
    stRefSelected->Draw();
    h1DataSelected->Draw("samee1");
    stRefSelected->GetXaxis()->SetTitle(titleX);
    lg_all->Draw();
    texSelected->Draw();
    texLumi->Draw();
    Utilities::saveCanvas(c1, "results/" + saveName  + "_lin_" + selected_name + "_" + det);

    // draw antiselected
    idealMaximumData = h1DataAntiselected->GetMaximum() + sqrt(h1DataAntiselected->GetMaximum())*2;
    if(stRefAntiselected->GetMaximum() < idealMaximumData ) stRefAntiselected->SetMaximum(idealMaximumData);
    stRefAntiselected->Draw();
    h1DataAntiselected->Draw("samee1");
    stRefAntiselected->GetXaxis()->SetTitle(titleX);
    lg_all->Draw();
    texAntiselected->Draw();
    texLumi->Draw();
    Utilities::saveCanvas(c1, "results/" + saveName  + "_lin_" + antiselected_name + "_" + det);

    delete texLumi;
    delete texSelected;
    delete texAntiselected;
    delete h1DataSelected;
    delete h1DataAntiselected;
    delete c1;
    delete stRefSelected;
    delete stRefAntiselected;


}

void plotResultsW(TString det, TString fileStamp, TString refFileStamp, const float &luminorm)
{

    gROOT->ProcessLine(".L tdrStyle.C");
    gROOT->ProcessLine("setTDRStyle()");

    // luminorm for 0.1pb-1 (is set to 1fb in the looper)
    //float luminorm = 0.0001;
    std::vector<DataSource> dataSources;
    std::vector<DataSource> refSources;

    dataSources.push_back (fH_WHUNT() );
    refSources.push_back( fH_WJETS() );
    refSources.push_back( fH_QCD30() );
    //refSources.push_back( fH_MINBIAS() );

    HistogramUtilities datah1(fileStamp + ".root", dataSources);
    HistogramUtilities refh1(refFileStamp + ".root", refSources, luminorm);

    // W studies related
    //

    // last argument of plotStack is rebin factor, cut val EB and cut val EE
    // in the case of topology "all" it uses cutValEE.

    //
    // Electrons
    //
    plotDataRefOverlayStack(datah1, refh1, "ele_selected_pt", "Electron p_{T} (GeV)", fileStamp, det, 8);
    plotDataRefOverlayStack(datah1, refh1, "ele_selected_pfmet", "Electron pfMET (GeV)", fileStamp, det, 4);
    plotDataRefOverlayStack(datah1, refh1, "ele_selected_tcmet", "Electron tcMET (GeV)", fileStamp, det, 4);
    plotDataRefOverlayStack(datah1, refh1, "ele_selected_ppfmet", "Electron Projected pfMET (GeV)", fileStamp, det, 4);
    plotDataRefOverlayStack(datah1, refh1, "ele_selected_ptcmet", "Electron Projected tcMET (GeV)", fileStamp, det, 4);
    plotDataRefOverlayStack(datah1, refh1, "ele_selected_tctransmass", "Transverse Mass w/tcMET (GeV)", fileStamp, det, 16);
    plotDataRefOverlayStack(datah1, refh1, "ele_selected_pftransmass", "Transverse Mass w/pfMET (GeV)", fileStamp, det, 16);
    plotDataRefOverlayStack(datah1, refh1, "ele_selected_d0corr", "d0(BS) (cm)", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "ele_selected_nmhits", "Number of Missing Hits", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "ele_selected_tcmetdphi", "dPhi(tcMET, Electron) (GeV)", fileStamp, det, 1);

    // N-1 distributions
    plotDataRefOverlayStack(datah1, refh1, "ele_nm1_tcmet", "NM1 tcMET (GeV)", fileStamp, det, 4, 20.0, 20.0);
    plotDataRefOverlayStack(datah1, refh1, "ele_nm1_pfmet", "NM1 pfMET (GeV)", fileStamp, det, 4, 20.0, 20.0);
    plotDataRefOverlayStack(datah1, refh1, "ele_nm1_iso", "NM1 RelIso (GeV)", fileStamp, det, 1, 0.10, 0.10);
    plotDataRefOverlayStack(datah1, refh1, "ele_nm1_jetveto", "NM1 Leading JPT p_{T} (GeV)", fileStamp, det, 1, 30.0, 30.0);
    plotDataRefOverlayStack(datah1, refh1, "ele_nm1_secondpt", "NM1 Second Electron p_{T} (GeV)", fileStamp, det, 1, 20.0, 20.0);
    plotDataRefOverlayStack(datah1, refh1, "ele_nm1_r19", "NM1 eMax/e5x5", fileStamp, det, 1, 0.95, 0.95);

    // N-1 distributions without data points
    plotStack(refh1, "ele_nm1_tcmet", "NM1 tcMET (GeV)", fileStamp + "_refonly", det, 4, 20.0, 20.0);
    plotStack(refh1, "ele_nm1_pfmet", "NM1 pfMET (GeV)", fileStamp + "_refonly", det, 4, 20.0, 20.0);
    plotStack(refh1, "ele_nm1_iso", "NM1 RelIso (GeV)", fileStamp + "_refonly", det, 1, 0.10, 0.10);
    plotStack(refh1, "ele_nm1_secondpt", "NM1 Second Electron p_{T} (GeV)", fileStamp + "_refonly", det, 1, 20.0, 20.0);
    plotStack(refh1, "ele_nm1_jetveto", "NM1 Leading JPT p_{T} (GeV)", fileStamp + "_refonly", det, 1, 30.0, 30.0);
    plotStack(refh1, "ele_nm1_r19", "NM1 eMax/e5x5", fileStamp + "_refonly", det, 1, 0.95, 0.95);

    // selection and anti-selection for electron id variable studies
    // barrel
    plotElectronIDStack(datah1, refh1, "sigmaIEtaIEta", "sigmaIEtaIEta", fileStamp, "eb", 2);
    plotElectronIDStack(datah1, refh1, "eOverPIn", "E/p_{IN}", fileStamp, "eb", 2);
    plotElectronIDStack(datah1, refh1, "hOverE", "H/E", fileStamp, "eb", 2);
    plotElectronIDStack(datah1, refh1, "e2x5MaxOver5x5", "E2x5Max/E5x5", fileStamp, "eb", 2);
    plotElectronIDStack(datah1, refh1, "dEtaIn", "dEtaIn", fileStamp, "eb", 2);
    plotElectronIDStack(datah1, refh1, "dPhiIn", "dPhiIn", fileStamp, "eb", 2);
    // endcap
    plotElectronIDStack(datah1, refh1, "sigmaIEtaIEta", "sigmaIEtaIEta", fileStamp, "ee", 2);
    plotElectronIDStack(datah1, refh1, "eOverPIn", "E/p_{IN}", fileStamp, "ee", 2);
    plotElectronIDStack(datah1, refh1, "hOverE", "H/E", fileStamp, "ee", 2);
    plotElectronIDStack(datah1, refh1, "e2x5MaxOver5x5", "E2x5Max/E5x5", fileStamp, "ee", 2);
    plotElectronIDStack(datah1, refh1, "dEtaIn", "dEtaIn", fileStamp, "ee", 2);
    plotElectronIDStack(datah1, refh1, "dPhiIn", "dPhiIn", fileStamp, "ee", 2);

    //
    // Muons
    //
    plotDataRefOverlayStack(datah1, refh1, "mu_selected_pt", "Muon p_{T} (GeV)", fileStamp, det, 8);
    plotDataRefOverlayStack(datah1, refh1, "mu_selected_pfmet", "Muon pfMET (GeV)", fileStamp, det, 4);
    plotDataRefOverlayStack(datah1, refh1, "mu_selected_tcmet", "Muon tcMET (GeV)", fileStamp, det, 4);
    plotDataRefOverlayStack(datah1, refh1, "mu_selected_ppfmet", "Muon Projected pfMET (GeV)", fileStamp, det, 4);
    plotDataRefOverlayStack(datah1, refh1, "mu_selected_ptcmet", "Muon Projected tcMET (GeV)", fileStamp, det, 4);
    plotDataRefOverlayStack(datah1, refh1, "mu_selected_tctransmass", "Transverse Mass w/tcMET (GeV)", fileStamp, det, 16);
    plotDataRefOverlayStack(datah1, refh1, "mu_selected_pftransmass", "Transverse Mass w/pfMET (GeV)", fileStamp, det, 16);
    plotDataRefOverlayStack(datah1, refh1, "mu_selected_d0corr", "d0(BS) (cm))", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "mu_selected_tcmetdphi", "dPhi(tcMET, Muon) (GeV)", fileStamp, det, 1);

    // N-1 distributions
    plotDataRefOverlayStack(datah1, refh1, "mu_nm1_tcmet", "NM1 tcMET (GeV)", fileStamp, det, 4, 20.0, 20.0);
    plotDataRefOverlayStack(datah1, refh1, "mu_nm1_pfmet", "NM1 pfMET (GeV)", fileStamp, det, 4, 20.0, 20.0);
    plotDataRefOverlayStack(datah1, refh1, "mu_nm1_iso", "NM1 RelIso (GeV)", fileStamp, det, 1, 0.10, 0.10);
    plotDataRefOverlayStack(datah1, refh1, "mu_nm1_secondpt", "NM1 Second Muon p_{T} (GeV)", fileStamp, det, 1, 20.0, 20.0);

    // N-1 distributions without data points
    plotStack(refh1, "mu_nm1_tcmet", "NM1 tcMET (GeV)", fileStamp + "_refonly", det, 4, 20.0, 20.0);
    plotStack(refh1, "mu_nm1_pfmet", "NM1 pfMET (GeV)", fileStamp + "_refonly", det, 4, 20.0, 20.0);
    plotStack(refh1, "mu_nm1_iso", "NM1 RelIso (GeV)", fileStamp + "_refonly", det, 1, 0.10, 0.10);
    plotStack(refh1, "mu_nm1_secondpt", "NM1 Second Muon p_{T} (GeV)", fileStamp + "_refonly", det, 1, 20.0, 20.0);





}



#include "plotResults.h"
#include "../../Tools/HistogramUtilities.h"
#include "../../Tools/AllDataSources.h"
#include "../../Tools/Utilities.h"
#include "../../Tools/histtools.cc"

#include "TROOT.h"

#include "TGraphAsymmErrors.h"
#include "TArrow.h"

// for 2_2_X
const static sources_t theMCSignal =  (1<<H_WJETS);

const static sources_t theMCBackground =  (1<<H_QCD30);

const static sources_t theMCSources = theMCSignal | theMCBackground;

const static sources_t theSources = theMCSignal | theMCBackground;

const static sources_t theDataSources = (1<<H_WHUNT);



TArrow *getArrow(THStack *st, TString det, float cutValEB, float cutValEE)
{

    float cutVal; 
    if (det == "eb") cutVal = cutValEB;
    else cutVal = cutValEE;
    TArrow *arr_cut = new TArrow(cutVal, st->GetMaximum()/2.0, cutVal, 0, 0.05, "|>");
    arr_cut->SetLineWidth(2.0);
    return arr_cut;
}


void plotStack(HistogramUtilities &h1, TString name, TString titleX, TString saveName, TString det, int rebin, float cutValEB, float cutValEE)
{

    THStack *st = h1.getStack(theSources, name, "", det, rebin);
    TLegend *lg_all = h1.getLegend(theSources, name, "", det);
    lg_all->SetX1(0.6);
    lg_all->SetX2(0.7);
    lg_all->SetY1(0.7);
    lg_all->SetY2(0.9);

    lg_all->SetTextSize(0.04);

    TCanvas *c1 = new TCanvas();
    c1->cd();
    st->Draw();
    st->GetXaxis()->SetTitle(titleX); 
    lg_all->Draw();
    if (det == "ee" && cutValEE != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    if (det == "eb" && cutValEB != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    Utilities::saveCanvas(c1, "results/" + saveName  + "_lin_" + name + "_" + det);

    c1->SetLogy();
    st->SetMinimum(1.0);
    st->Draw();	
    lg_all->Draw();
    if (det == "ee" && cutValEE != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    if (det == "eb" && cutValEB != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    Utilities::saveCanvas(c1, "results/" + saveName  + "_log_" + name + "_" + det);

    delete c1;
    delete st;
    delete lg_all;

}


void plotDataRefOverlayStack(HistogramUtilities &hData, HistogramUtilities &hRef, 
        TString name, TString titleX, TString saveName, TString det, int rebin)
{

    THStack *stRef = hRef.getStack(theMCSources, name, "", det, rebin);
    TLegend *lg_all = hRef.getLegend(theMCSources, name, "", det);
    lg_all->SetX1(0.6);
    lg_all->SetX2(0.7);
    lg_all->SetY1(0.7);
    lg_all->SetY2(0.9);
    lg_all->SetTextSize(0.04);

    TH1F *h1Data =  hData.getHistogram(theDataSources, name, "", det, rebin);
    h1Data->GetXaxis()->SetTitle("");
    h1Data->GetYaxis()->SetTitle("");
    h1Data->GetXaxis()->SetLabelSize(0);
    h1Data->GetYaxis()->SetLabelSize(0);
    h1Data->SetMarkerStyle(20);
    h1Data->SetMarkerSize(1.5);
    h1Data->SetLineWidth(1.5);
    h1Data->SetMarkerColor(kRed);
    h1Data->SetLineColor(kRed);
    lg_all->AddEntry(h1Data,"Data", "lep"); 
    if(stRef->GetMaximum() < h1Data->GetMaximum() )
        stRef->SetMaximum(h1Data->GetMaximum()*2.0);

    TCanvas *c1 = new TCanvas();
    c1->cd();
    stRef->Draw();
    h1Data->Draw("samee1");
    stRef->GetXaxis()->SetTitle(titleX); 
    lg_all->Draw();
    Utilities::saveCanvas(c1, "results/" + saveName  + "_lin_" + name + "_" + det);

    c1->SetLogy();
    stRef->SetMinimum(1.0);
    stRef->Draw();	
    h1Data->Draw("samee1");
    lg_all->Draw();
    Utilities::saveCanvas(c1, "results/" + saveName  + "_log_" + name + "_" + det);

    delete c1;
    delete stRef;
    delete lg_all;

}


void plotResultsW(TString det, TString fileStamp, TString refFileStamp, float & luminorm)
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

    // last argument of plotStack is rebin factor

    // Electrons
    plotDataRefOverlayStack(datah1, refh1, "ele_selected_pt", "Electron p_{T} (GeV)", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "ele_selected_ppfmet", "Electron Projected pfMET (GeV)", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "ele_selected_ptcmet", "Electron Projected tcMET (GeV)", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "ele_selected_tctransmass", "Transverse Mass w/tcMET (GeV)", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "ele_selected_pftransmass", "Transverse Mass w/pfMET (GeV)", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "ele_selected_d0corr", "d0(BS) (cm)", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "ele_selected_nmhits", "Number of Missing Hits", fileStamp, det, 1);

    // N-1 distributions
    plotDataRefOverlayStack(datah1, refh1, "ele_nm1_tcmet", "NM1 tcMET (GeV)", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "ele_nm1_pfmet", "NM1 pfMET (GeV)", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "ele_nm1_iso", "NM1 RelIso (GeV)", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "ele_nm1_jetveto", "NM1 Leading JPT p_{T} (GeV)", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "ele_nm1_secondpt", "NM1 Second Electron p_{T} (GeV)", fileStamp, det, 1);


    // Muons
    plotDataRefOverlayStack(datah1, refh1, "mu_selected_pt", "Muon p_{T} (GeV)", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "mu_selected_ppfmet", "Muon Projected pfMET (GeV)", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "mu_selected_ptcmet", "Muon Projected tcMET (GeV)", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "mu_selected_tctransmass", "Transverse Mass w/tcMET (GeV)", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "mu_selected_pftransmass", "Transverse Mass w/pfMET (GeV)", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "mu_selected_d0corr", "d0(BS) (cm))", fileStamp, det, 1);

    // N-1 distributions
    plotDataRefOverlayStack(datah1, refh1, "mu_nm1_tcmet", "NM1 tcMET (GeV)", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "mu_nm1_pfmet", "NM1 pfMET (GeV)", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "mu_nm1_iso", "NM1 RelIso (GeV)", fileStamp, det, 1);
    plotDataRefOverlayStack(datah1, refh1, "mu_nm1_secondpt", "NM1 Second Muon p_{T} (GeV)", fileStamp, det, 1);


}



#include "plotResults.h"
#include "../../../Tools/HistogramUtilities.h"
#include "../../../Tools/AllDataSources.h"
#include "../../../Tools/Utilities.h"
#include "../../../Tools/histtools.cc"

#include "TROOT.h"

#include "TGraphAsymmErrors.h"
#include "TArrow.h"
#include "TLatex.h"

const static sources_t theMCSignal =  (1<<H_TTBAR);
const static sources_t theMCBackground = (1<<H_ZJETS) | (1<<H_WJETS) | (1<<H_QCD30);
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
    st->Draw("hist");
    st->GetXaxis()->SetTitle(titleX); 
    lg_all->Draw();
    if (det == "ee" && cutValEE != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    if (det == "eb" && cutValEB != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    if (det == "all" && cutValEE != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    Utilities::saveCanvas(c1, "../results/" + saveName  + "_lin_" + name + "_" + det);

    c1->SetLogy();
    st->SetMinimum(0.10);
    st->Draw("hist");	
    lg_all->Draw();
    if (det == "ee" && cutValEE != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    if (det == "eb" && cutValEB != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    if (det == "all" && cutValEE != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    Utilities::saveCanvas(c1, "../results/" + saveName  + "_log_" + name + "_" + det);

    delete c1;
    delete st;
    delete lg_all;

}


void plotDataRefOverlayStack(HistogramUtilities &hData, HistogramUtilities &hRef, 
        TString name, TString titleX, TString saveName, TString det, TString norm, int rebin,  float cutValEB, float cutValEE)
{

    TLatex *   tex = new TLatex(0.55,0.83,"MC Norm (" + norm + ")");
    tex->SetNDC();
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    TLatex *   tex2 = new TLatex(0.55,0.88,"CMS Preliminary 2010");
    tex2->SetNDC();
    tex2->SetTextSize(0.04);
    tex2->SetLineWidth(2);
    TLatex *   tex3 = new TLatex(0.25, 0.88,"Final State: " + det);
    tex3->SetNDC();
    tex3->SetTextSize(0.04);
    tex3->SetLineWidth(2);

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
    h1Data->SetMarkerColor(kBlack);
    h1Data->SetLineColor(kBlack);
    lg_all->AddEntry(h1Data,"Data", "lep"); 

    // ideal y range should be at least 2 sigma greater than highest data point
    float idealMaximumData = h1Data->GetMaximum() + sqrt(h1Data->GetMaximum())*2;
    if(stRef->GetMaximum() < idealMaximumData ) stRef->SetMaximum(idealMaximumData);

    TCanvas *c1 = new TCanvas();
    c1->cd();
    stRef->Draw("HIST");
    h1Data->Draw("samee1");
    stRef->GetXaxis()->SetTitle(titleX); 
    lg_all->Draw();
    tex->Draw();
    tex2->Draw();
    tex3->Draw();
    if (det == "ee" && cutValEE != -1) getArrow(stRef, det, cutValEB, cutValEE, idealMaximumData/1.2)->Draw();
    if (det == "eb" && cutValEB != -1) getArrow(stRef, det, cutValEB, cutValEE, idealMaximumData/1.2)->Draw();
    if (det == "all" && cutValEE != -1) getArrow(stRef, det, cutValEB, cutValEE, idealMaximumData/1.2)->Draw();
    Utilities::saveCanvas(c1, "../results/" + saveName  + "_lin_" + name + "_" + det);

    c1->SetLogy();
    stRef->SetMinimum(0.10);
    stRef->Draw("HIST");	
    h1Data->Draw("samee1");
    lg_all->Draw();
    tex->Draw();
    tex2->Draw();
    tex3->Draw();
    if (det == "ee" && cutValEE != -1) getArrow(stRef, det, cutValEB, cutValEE, idealMaximumData/1.2)->Draw();
    if (det == "eb" && cutValEB != -1) getArrow(stRef, det, cutValEB, cutValEE, idealMaximumData/1.2)->Draw();
    if (det == "all" && cutValEE != -1) getArrow(stRef, det, cutValEB, cutValEE, idealMaximumData/1.2)->Draw();
    Utilities::saveCanvas(c1, "../results/" + saveName  + "_log_" + name + "_" + det);

    delete tex;
    delete tex2;
    delete c1;
    delete stRef;
    delete lg_all;

}

void plotResults(TString det, TString fileStamp, TString refFileStamp, TString norm, const float &luminorm)
{

    gROOT->ProcessLine(".L ../tdrStyle.C");
    gROOT->ProcessLine("setTDRStyle()");

    // luminorm for 0.1pb-1 (is set to 1fb in the looper)
    //float luminorm = 0.0001;
    std::vector<DataSource> dataSources;
    std::vector<DataSource> refSources;
    dataSources.push_back (fH_WHUNT() );
    refSources.push_back( fH_TTBAR() );
    refSources.push_back( fH_ZJETS() );
    refSources.push_back( fH_QCD30() );
    refSources.push_back( fH_WJETS() );

    HistogramUtilities datah1("../" + fileStamp + ".root", dataSources);
    HistogramUtilities refh1("../" + refFileStamp + ".root", refSources, luminorm);

    //
    // Electrons
    //

    // before any selections
    plotDataRefOverlayStack(datah1, refh1, "hyp_tcmet_allj", "tcMET (GeV)", fileStamp, det, norm, 1);
    plotDataRefOverlayStack(datah1, refh1, "hyp_ll_pt_allj", "pT (LL) (GeV)", fileStamp, det, norm, 2);
    plotDataRefOverlayStack(datah1, refh1, "hyp_lt_pt_allj", "pT (LT) (GeV)", fileStamp, det, norm, 2);
    plotDataRefOverlayStack(datah1, refh1, "hyp_ll_d0corr_allj", "d0corr (LL) cm", fileStamp, det, norm, 2);
    plotDataRefOverlayStack(datah1, refh1, "hyp_lt_d0corr_allj", "d0corr (LT) cm", fileStamp, det, norm, 2);
    plotDataRefOverlayStack(datah1, refh1, "hyp_allmuon_d0corr_allj", "d0corr (all muons) cm", fileStamp, det, norm, 4);
    plotDataRefOverlayStack(datah1, refh1, "hyp_glbtrkmuon_d0corr_allj", "d0corr (glb+trk muons) cm", fileStamp, det, norm, 4);



    // preselection
    plotDataRefOverlayStack(datah1, refh1, "hyp_presel_njets", "Number of Jets", fileStamp, det, norm, 1);
    plotDataRefOverlayStack(datah1, refh1, "hyp_presel_tcmet_allj", "tcMET (GeV)", fileStamp, det, norm, 1);
    plotDataRefOverlayStack(datah1, refh1, "hyp_presel_mll_allj", "M_{ll} (GeV)", fileStamp, det, norm, 1);
    plotDataRefOverlayStack(datah1, refh1, "hyp_presel_pt_allj", "pT_{ll} (GeV)", fileStamp, det, norm, 4);
    plotDataRefOverlayStack(datah1, refh1, "hyp_presel_lt_iso_allj", "RelIso (LT) (GeV)", fileStamp, det, norm, 2);
    plotDataRefOverlayStack(datah1, refh1, "hyp_presel_ll_iso_allj", "RelIso (LL) (GeV)", fileStamp, det, norm, 2);
    plotDataRefOverlayStack(datah1, refh1, "hyp_presel_ll_pt_allj", "pT (LL) (GeV)", fileStamp, det, norm, 2);
    plotDataRefOverlayStack(datah1, refh1, "hyp_presel_lt_pt_allj", "pT (LT) (GeV)", fileStamp, det, norm, 2);
    plotDataRefOverlayStack(datah1, refh1, "hyp_presel_ll_d0corr_allj", "d0corr (LL) cm", fileStamp, det, norm, 2);
    plotDataRefOverlayStack(datah1, refh1, "hyp_presel_lt_d0corr_allj", "d0corr (LT) cm", fileStamp, det, norm, 2);

    // preselection + loose iso
    plotDataRefOverlayStack(datah1, refh1, "hyp_presel_looseiso_njets", "Number of Jets", fileStamp, det, norm, 1);
    plotDataRefOverlayStack(datah1, refh1, "hyp_presel_looseiso_tcmet_allj", "tcMET (GeV)", fileStamp, det, norm, 1);
    plotDataRefOverlayStack(datah1, refh1, "hyp_presel_looseiso_mll_allj", "M_{ll} (GeV)", fileStamp, det, norm, 1);
    plotDataRefOverlayStack(datah1, refh1, "hyp_presel_looseiso_pt_allj", "pT_{ll} (GeV)", fileStamp, det, norm, 4);


}


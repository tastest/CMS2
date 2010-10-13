
#include "plotResults.h"
#include "../Tools/HistogramUtilities.h"
#include "../Tools/AllDataSources.h"
#include "../Tools/Utilities.h"
#include "../Tools/histtools.cc"

#include "TROOT.h"

#include "TGraphAsymmErrors.h"
#include "TArrow.h"
#include "TLatex.h"

const static sources_t theMCSignal =  (1ll<<H_TTBAR);
const static sources_t theMCBackground = (1ll<<H_DYEE) | (1ll<<H_DYMM) | (1ll<<H_DYTT) |  (1ll<<H_WJETS);
const static sources_t theMCSources = theMCSignal | theMCBackground;
const static sources_t theDataSources = (1ll<<H_DATA);

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


void plotStack(HistogramUtilities &h1, TString name, TString titleX, TString saveName, TString det, TString norm, int rebin, float cutValEB, float cutValEE)
{

    TLatex *   tex = new TLatex(0.55,0.83,"MC Norm (" + norm + ")");
    tex->SetNDC();
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    TLatex *   tex2 = new TLatex(0.55,0.88,"CMS Preliminary 2010");
    tex2->SetNDC();
    tex2->SetTextSize(0.04);
    tex2->SetLineWidth(2);

    THStack *st = h1.getStack(theMCSources, name, "", det, rebin);
    TLegend *lg_all = h1.getLegend(theMCSources, name, "", det);
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
    tex->Draw();
    tex2->Draw();
    if (det == "ee" && cutValEE != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    if (det == "eb" && cutValEB != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    if (det == "all" && cutValEE != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    Utilities::saveCanvas(c1, "results/" + saveName  + "_lin_" + name + "_" + det);

    c1->SetLogy();
    st->SetMinimum(0.10);
    st->Draw("hist");	
    lg_all->Draw();
    tex->Draw();
    tex2->Draw();
    if (det == "ee" && cutValEE != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    if (det == "eb" && cutValEB != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    if (det == "all" && cutValEE != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
    Utilities::saveCanvas(c1, "results/" + saveName  + "_log_" + name + "_" + det);

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
    stRef->Draw("HIST");
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
    stRef->Draw("HIST");	
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

void plotResults(TString hyp, TString version, TString fileStamp, TString refFileStamp, TString norm, const float &luminorm)
{

    gROOT->ProcessLine(".L tdrStyle.C");
    gROOT->ProcessLine("setTDRStyle()");

    // luminorm for 0.1pb-1 (is set to 1fb in the looper)
    //float luminorm = 0.0001;
    std::vector<DataSource> dataSources;
    std::vector<DataSource> refSources;

    dataSources.push_back (fH_DATA() );
    refSources.push_back( fH_DYEE() );
    refSources.push_back( fH_DYMM() );
    refSources.push_back( fH_DYTT() );
    refSources.push_back( fH_TTBAR() );
    refSources.push_back( fH_WJETS() );

    //HistogramUtilities datah1(fileStamp + ".root", dataSources);
    HistogramUtilities refh1(refFileStamp + ".root", refSources, luminorm);

    // last argument of plotStack is rebin factor, cut val EB and cut val EE
    // in the case of topology "all" it uses cutValEE.

    plotStack(refh1, "hnSumPTall", "M_{Eff} (GeV)", version + "_mc", hyp, norm, 1);

}

void makeTables(TString fileStamp, TString refFileStamp, TString norm, const float &luminorm)
{

    std::vector<DataSource> dataSources;
    std::vector<DataSource> refSources;

    dataSources.push_back (fH_DATA() );
    refSources.push_back( fH_DYEE() );
    refSources.push_back( fH_DYMM() );
    refSources.push_back( fH_DYTT() );
    refSources.push_back( fH_TTBAR() );
    refSources.push_back( fH_WJETS() );

    //HistogramUtilities datah1(fileStamp + ".root", dataSources);
    HistogramUtilities refh1(refFileStamp + ".root", refSources, luminorm);

    static const char dilepton_hypo_names[][128] = { "all", "mm", "em", "ee" };

    // print table
    std::cout << "\\begin{table}[ht]" << std::endl;
    std::cout << "\\caption{Yields for $" <<  norm << "$}" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|l*{5}{|c}|r|}\\hline" << std::endl;
    std::cout << "Source & "    << dilepton_hypo_names[0] << " & " << dilepton_hypo_names[1] << " & " 
                                << dilepton_hypo_names[2] << " & " << dilepton_hypo_names[3] << "\\\\ \\hline" << std::endl;

    for (size_t i = 0; i < refSources.size(); ++i) {
        std::cout << refSources[i].getName();
        for (unsigned int h = 0; h < 4; ++h) {
            TH1F* h_temp = refh1.getHistogram((1ll<<refSources[i].getSource()), "hnSumPTall", "", dilepton_hypo_names[h]);
            // note that overflow bin is already moved to last bit by HistogramUtilities::getHistogram
            std::cout << " & " << h_temp->Integral(0, h_temp->GetNbinsX());
            delete h_temp;
        }
        std::cout << " \\\\ \\hline" << std::endl;
    }    

    std::cout <<"\\end{tabular}" << std::endl;
    std::cout <<"\\end{center}" << std::endl;
    std::cout << "\\end{table}" << std::endl;


}



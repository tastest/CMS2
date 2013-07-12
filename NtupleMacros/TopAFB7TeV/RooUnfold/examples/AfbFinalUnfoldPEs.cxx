#include <iostream>
#include <fstream>
#include "AfbFinalUnfold.h"

#include "TROOT.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TMatrixD.h"
#include "TChain.h"
#include "TLegend.h"
#include "TColor.h"
#include "THStack.h"
#include "TCut.h"
#include "TPaveText.h"
#include "TLatex.h"

#include "TUnfold.h"

#include "src/RooUnfoldResponse.h"
#include "src/RooUnfoldBayes.h"
#include "src/RooUnfoldSvd.h"
#include "src/RooUnfoldTUnfold.h"

#include "tdrstyle.C"

using std::cout;
using std::endl;


//==============================================================================
// Global definitions
//==============================================================================

// 0=SVD, 1=TUnfold via RooUnfold, 2=TUnfold
int unfoldingType = 0;

TString Region = "";
Int_t kterm = 3;
Double_t tau = 1E-4;
Int_t nVars = 8;
Int_t includeSys = 0;


void AfbUnfoldExample(double scalettdil = 1., double scalettotr = 1., double scalewjets = 1., double scaleDY = 1., double scaletw = 1., double scaleVV = 1. )
{

    setTDRStyle();
    gStyle->SetOptFit();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    cout.precision(3);

    TString summary_name = "summary_1Dunfolding";

    if (!(scalettotr == 1. && scalewjets == 1. && scaleDY == 1. && scaletw == 1. && scaleVV == 1.))  summary_name = Form("summary_1Dunfolding_%i_%i_%i_%i_%i", int(10.*scalettotr + 0.5), int(10.*scalewjets + 0.5), int(10.*scaleDY + 0.5), int(10.*scaletw + 0.5), int(10.*scaleVV + 0.5));

    ofstream myfile;
    myfile.open (summary_name + ".txt");
    cout.rdbuf(myfile.rdbuf());

    // OGU 130516: add second output txt file with format easier to be pasted into google docs
    ofstream second_output_file;
    second_output_file.open(summary_name + "_formated.txt");

    TRandom3 *random = new TRandom3();
    random->SetSeed(5);


    const int nBkg = 7;
    TString path = "../";
    TString bkgroot[nBkg] = {"ttotr.root", "wjets.root", "DYee.root", "DYmm.root", "DYtautau.root", "tw.root", "VV.root"};
    double bkgSF[nBkg] = {scalettotr, scalewjets, scaleDY, scaleDY, scaleDY, scaletw, scaleVV};

    Float_t observable, observable_gen, ttmass, ttRapidity, tmass;
    Float_t observableMinus, observableMinus_gen;
    Double_t weight;
    Int_t evt;
    Int_t Nsolns;
    const int NPEs = 900;
    Int_t PEsize = 1000;

    for (Int_t iVar = 0; iVar < nVars; iVar++)
    {

        Initialize1DBinning(iVar);
        bool combineLepMinus = acceptanceName == "lepCosTheta" ? true : false;

        TH1D *hData = new TH1D ("Data_BkgSub", "Data with background subtracted",    nbins1D, xbins1D);
        TH1D *hBkg = new TH1D ("Background",  "Background",    nbins1D, xbins1D);
        TH1D *hData_unfolded = new TH1D ("Data_Unfold", "Data with background subtracted and unfolded", nbins1D, xbins1D);

        TH1D *hTrue = new TH1D ("true", "Truth",    nbins1D, xbins1D);
        TH1D *hMeas = new TH1D ("meas", "Measured", nbins1D, xbins1D);

        TH2D *hTrue_vs_Meas = new TH2D ("true_vs_meas", "True vs Measured", nbins1D, xbins1D, nbins1D, xbins1D);

        TH1D *hData_bkgSub;

        TH1D *hPseudoData[NPEs];
        TH1D *hPseudoData_unfolded[NPEs];
        TMatrixD m_unfoldE_PEs[NPEs];
        TMatrixD m_correctE_PEs[NPEs];

        for (int iPE = 0; iPE < NPEs; ++iPE)
        {
            hPseudoData[iPE] = new TH1D ("", "", nbins1D, xbins1D);
            hPseudoData_unfolded[iPE] = new TH1D ("", "", nbins1D, xbins1D);
            m_unfoldE_PEs[iPE] = new TMatrixD(nbins1D, nbins1D);
            m_correctE_PEs[iPE] = new TMatrixD(nbins1D, nbins1D);
        }

        hData->Sumw2();
        hBkg->Sumw2();
        hData_unfolded->Sumw2();
        hTrue->Sumw2();
        hMeas->Sumw2();
        hTrue_vs_Meas->Sumw2();

        TMatrixD m_unfoldE (nbins1D, nbins1D);
        TMatrixD m_correctE(nbins1D, nbins1D);


        //  Now test with data and with BKG subtraction

        TChain *ch_bkg[nBkg];
        TChain *ch_top = new TChain("tree");

        TChain *ch_data = new TChain("tree");


        ch_data->Add(path + "data.root");

        ch_top->Add(path + "ttdil.root");

        for (int iBkg = 0; iBkg < nBkg; ++iBkg)
        {
            ch_bkg[iBkg] = new TChain("tree");
            ch_bkg[iBkg]->Add(path + bkgroot[iBkg]);
        }


        ch_data->SetBranchAddress(observablename,    &observable);
        if ( combineLepMinus ) ch_data->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
        ch_data->SetBranchAddress("weight", &weight);
        ch_data->SetBranchAddress("Nsolns", &Nsolns);
        ch_data->SetBranchAddress("tt_mass", &ttmass);
        ch_data->SetBranchAddress("ttRapidity", &ttRapidity);
        ch_data->SetBranchAddress("t_mass", &tmass);


        for (Int_t i = 0; i < ch_data->GetEntries(); i++)
        {
            ch_data->GetEntry(i);
            if ( (acceptanceName == "lepChargeAsym") || (acceptanceName == "lepAzimAsym") || (acceptanceName == "lepAzimAsym2") )
            {
                // leptonic asymmetries don't need valid top mass solution
                fillUnderOverFlow(hData, observable, weight, Nsolns);
            }
            else
            {
                if ( ttmass > 0 )
                {
                    // asymmetries with top properties are required to have a valid top mass solution
                    fillUnderOverFlow(hData, observable, weight, Nsolns);
                }
            }
            if (combineLepMinus)
            {
                // combine plus and minus
                if ( (acceptanceName == "lepChargeAsym") || (acceptanceName == "lepAzimAsym") || (acceptanceName == "lepAzimAsym2") )
                {
                    // leptonic asymmetries don't need valid top mass solution
                    fillUnderOverFlow(hData, observableMinus, weight, Nsolns);
                }
                else
                {
                    if ( ttmass > 0 )
                    {
                        // asymmetries with top properties are required to have a valid top mass solution
                        fillUnderOverFlow(hData, observableMinus, weight, Nsolns);
                    }
                }
            }
        }


        for (int iBkg = 0; iBkg < nBkg; ++iBkg)
        {

            ch_bkg[iBkg]->SetBranchAddress(observablename,    &observable);
            if ( combineLepMinus ) ch_bkg[iBkg]->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
            ch_bkg[iBkg]->SetBranchAddress("weight", &weight);
            ch_bkg[iBkg]->SetBranchAddress("Nsolns", &Nsolns);
            ch_bkg[iBkg]->SetBranchAddress("tt_mass", &ttmass);
            ch_bkg[iBkg]->SetBranchAddress("ttRapidity", &ttRapidity);
            ch_bkg[iBkg]->SetBranchAddress("t_mass", &tmass);

            for (Int_t i = 0; i < ch_bkg[iBkg]->GetEntries(); i++)
            {
                ch_bkg[iBkg]->GetEntry(i);
                weight *= bkgSF[iBkg];
                if ( (acceptanceName == "lepChargeAsym") || (acceptanceName == "lepAzimAsym") || (acceptanceName == "lepAzimAsym2") )
                {
                    // leptonic asymmetries don't need valid top mass solution
                    fillUnderOverFlow(hBkg, observable, weight, Nsolns);
                }
                else
                {
                    if ( ttmass > 0 )
                    {
                        // asymmetries with top properties are required to have a valid top mass solution
                        fillUnderOverFlow(hBkg, observable, weight, Nsolns);
                    }
                }
                if (combineLepMinus)
                {
                    if ( (acceptanceName == "lepChargeAsym") || (acceptanceName == "lepAzimAsym") || (acceptanceName == "lepAzimAsym2") )
                    {
                        // leptonic asymmetries don't need valid top mass solution
                        fillUnderOverFlow(hBkg, observableMinus, weight, Nsolns);
                    }
                    else
                    {
                        if ( ttmass > 0 )
                        {
                            // asymmetries with top properties are required to have a valid top mass solution
                            fillUnderOverFlow(hBkg, observableMinus, weight, Nsolns);
                        }
                    }
                }
            }

        }

        ch_top->SetBranchAddress(observablename,    &observable);
        ch_top->SetBranchAddress(observablename + "_gen", &observable_gen);
        if (observablename == "lep_azimuthal_asymmetry2") ch_top->SetBranchAddress("lep_azimuthal_asymmetry_gen2", &observable_gen);
        if ( combineLepMinus ) ch_top->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
        if ( combineLepMinus ) ch_top->SetBranchAddress("lepMinus_costheta_cms_gen",    &observableMinus_gen);
        ch_top->SetBranchAddress("weight", &weight);
        ch_top->SetBranchAddress("Nsolns", &Nsolns);
        ch_top->SetBranchAddress("tt_mass", &ttmass);
        ch_top->SetBranchAddress("ttRapidity", &ttRapidity);
        ch_top->SetBranchAddress("t_mass", &tmass);
        ch_top->SetBranchAddress("evt", &evt);

        Int_t prevevt = -999;
        Int_t i_ev = -1;

        for (Int_t i = 0; i < ch_top->GetEntries(); i++)
        {
            ch_top->GetEntry(i);
            weight *= scalettdil;

            if (evt != prevevt) i_ev++;
            //if(evt!=prevevt) cout<<i_ev<<" "<<int(i_ev/PEsize)<<" "<<Nsolns<<endl;


            if ( (acceptanceName == "lepChargeAsym") || (acceptanceName == "lepAzimAsym") || (acceptanceName == "lepAzimAsym2") )
            {
                fillUnderOverFlow(hMeas, observable, weight, Nsolns);
                fillUnderOverFlow(hTrue, observable_gen, weight, Nsolns);
                fillUnderOverFlow(hTrue_vs_Meas, observable, observable_gen, weight, Nsolns);
                if (int(i_ev / PEsize) < NPEs) fillUnderOverFlow(hPseudoData[int(i_ev / PEsize)], observable, weight, Nsolns);
                if ( combineLepMinus )
                {
                    fillUnderOverFlow(hMeas, observableMinus, weight, Nsolns);
                    fillUnderOverFlow(hTrue, observableMinus_gen, weight, Nsolns);
                    fillUnderOverFlow(hTrue_vs_Meas, observableMinus, observableMinus_gen, weight, Nsolns);
                    if (int(i_ev / PEsize) < NPEs) fillUnderOverFlow(hPseudoData[int(i_ev / PEsize)], observableMinus, weight, Nsolns);
                }
            }
            else
            {
                if ( ttmass > 0 )
                {
                    fillUnderOverFlow(hMeas, observable, weight, Nsolns);
                    fillUnderOverFlow(hTrue, observable_gen, weight, Nsolns);
                    fillUnderOverFlow(hTrue_vs_Meas, observable, observable_gen, weight, Nsolns);
                    if (int(i_ev / PEsize) < NPEs) fillUnderOverFlow(hPseudoData[int(i_ev / PEsize)], observable, weight, Nsolns);
                    if ( combineLepMinus )
                    {
                        fillUnderOverFlow(hMeas, observableMinus, weight, Nsolns);
                        fillUnderOverFlow(hTrue, observableMinus_gen, weight, Nsolns);
                        fillUnderOverFlow(hTrue_vs_Meas, observableMinus, observableMinus_gen, weight, Nsolns);
                        if (int(i_ev / PEsize) < NPEs) fillUnderOverFlow(hPseudoData[int(i_ev / PEsize)], observableMinus, weight, Nsolns);
                    }
                }
            }

            prevevt = evt;
        }

        Int_t lastPE = int(i_ev / PEsize);

        cout << "last PE: " << lastPE << " " << i_ev << endl;

        RooUnfoldResponse response (hMeas, hTrue, hTrue_vs_Meas);

        hData_bkgSub = (TH1D *) hData->Clone();
        hData_bkgSub->Add(hBkg, -1.0);

        if (unfoldingType == 0)
        {
            RooUnfoldSvd unfold_svd (&response, hData_bkgSub, kterm);
            unfold_svd.Setup(&response, hData_bkgSub);
            unfold_svd.IncludeSystematics(includeSys);
            hData_unfolded = (TH1D *) unfold_svd.Hreco();
            m_unfoldE = unfold_svd.Ereco();

            for (int iPE = 0; iPE < lastPE; ++iPE)
            {
                RooUnfoldSvd unfold_svd (&response, hPseudoData[iPE], kterm);
                unfold_svd.Setup(&response, hPseudoData[iPE]);
                unfold_svd.IncludeSystematics(includeSys);
                hPseudoData_unfolded[iPE] = (TH1D *) unfold_svd.Hreco();
                m_unfoldE_PEs[iPE] = unfold_svd.Ereco();
            }
        }
        else if (unfoldingType == 1)
        {
            RooUnfoldTUnfold unfold_rooTUnfold (&response, hData_bkgSub, TUnfold::kRegModeCurvature);
            unfold_rooTUnfold.Setup(&response, hData_bkgSub);
            unfold_rooTUnfold.FixTau(tau);
            unfold_rooTUnfold.IncludeSystematics(includeSys);
            hData_unfolded = (TH1D *) unfold_rooTUnfold.Hreco();
            m_unfoldE = unfold_rooTUnfold.Ereco();
        }
        else if (unfoldingType == 2)
        {
            TUnfold unfold_TUnfold (hTrue_vs_Meas, TUnfold::kHistMapOutputVert, TUnfold::kRegModeCurvature);
            unfold_TUnfold.SetInput(hMeas);
            //Double_t biasScale=1.0;
            unfold_TUnfold.SetBias(hTrue);
            unfold_TUnfold.DoUnfold(tau);
            unfold_TUnfold.GetOutput(hData_unfolded);


            TH2D *ematrix = unfold_TUnfold.GetEmatrix("ematrix", "error matrix", 0, 0);
            for (Int_t cmi = 0; cmi < nbins1D; cmi++)
            {
                for (Int_t cmj = 0; cmj < nbins1D; cmj++)
                {
                    m_unfoldE(cmi, cmj) = ematrix->GetBinContent(cmi + 1, cmj + 1);
                }
            }
        }
        else cout << "Unfolding TYPE not Specified" << "\n";


        if (unfoldingType == 0)
        {
            TCanvas *c_d = new TCanvas("c_d", "c_d", 500, 500);
            TH1D *dvec = unfold_svd.Impl()->GetD();
            dvec->Draw();
            c_d->SetLogy();
            c_d->SaveAs("D_" + acceptanceName + Region + ".pdf");
        }

        TCanvas *c_resp = new TCanvas("c_resp", "c_resp");
        TH2D *hResp = (TH2D *) response.Hresponse();
        gStyle->SetPalette(1);
        hResp->GetXaxis()->SetTitle(xaxislabel + "_{gen}");
        hResp->GetYaxis()->SetTitle(xaxislabel);
        hResp->Draw("COLZ");
        c_resp->SaveAs("Response_" + acceptanceName + Region + ".pdf");

        TFile *file = new TFile("../acceptance/mcnlo/accept_" + acceptanceName + ".root");
        TH1D *acceptM = (TH1D *) file->Get("accept_" + acceptanceName);
        acceptM->Scale(1.0 / acceptM->Integral());

        TH1D *denominatorM = (TH1D *) file->Get("denominator_" + acceptanceName);

        for (Int_t i = 1; i <= nbins1D; i++)
        {

            if (acceptM->GetBinContent(i) != 0)
            {
                hData_unfolded->SetBinContent(i, hData_unfolded->GetBinContent(i) * 1.0 / acceptM->GetBinContent(i));
                hData_unfolded->SetBinError  (i, hData_unfolded->GetBinError  (i) * 1.0 / acceptM->GetBinContent(i));
                for (int iPE = 0; iPE < lastPE; ++iPE)
                {
                    hPseudoData_unfolded[iPE]->SetBinContent(i, hPseudoData_unfolded[iPE]->GetBinContent(i) * 1.0 / acceptM->GetBinContent(i));
                    hPseudoData_unfolded[iPE]->SetBinError  (i, hPseudoData_unfolded[iPE]->GetBinError  (i) * 1.0 / acceptM->GetBinContent(i));
                }
            }

            if (acceptM->GetBinContent(i) != 0)
            {
                hTrue->SetBinContent(i, hTrue->GetBinContent(i) * 1.0 / acceptM->GetBinContent(i));
                hTrue->SetBinError  (i, hTrue->GetBinError(i)  * 1.0 / acceptM->GetBinContent(i));
            }
        }

        for (int l = 0; l < nbins1D; l++)
        {
            for (int j = 0; j < nbins1D; j++)
            {
                double corr = 1.0 / ( acceptM->GetBinContent(l + 1) * acceptM->GetBinContent(j + 1) );
                //corr = corr * pow(xsection / dataIntegral,2) ;
                m_correctE(l, j) = m_unfoldE(l, j) * corr;
            }
        }


        for (int iPE = 0; iPE < lastPE; ++iPE)
        {
            TMatrixD m_unfoldE_temp =  m_unfoldE_PEs[iPE];
            TMatrixD m_correctE_temp =  m_correctE_PEs[iPE];
            for (int l = 0; l < nbins1D; l++)
            {
                for (int j = 0; j < nbins1D; j++)
                {
                    double corr = 1.0 / ( acceptM->GetBinContent(l + 1) * acceptM->GetBinContent(j + 1) );
                    m_correctE_temp(l, j) = m_unfoldE_temp(l, j) * corr;
                }
            }
            m_correctE_PEs[iPE] = m_correctE_temp;
        }


        //==================================================================
        // ============== Print the asymmetry =============================
        cout << "========= Variable: " << acceptanceName << "===================\n";

        Float_t Afb, AfbErr;

        GetAfb(hData, Afb, AfbErr);
        cout << " Data: " << Afb << " +/-  " << AfbErr << "\n";

        GetAfb(hTrue, Afb, AfbErr);
        cout << " True Top: " << Afb << " +/-  " << AfbErr << "\n";

        GetCorrectedAfb(hData_unfolded, m_correctE, Afb, AfbErr);
        cout << " Unfolded: " << Afb << " +/-  " << AfbErr << "\n";
        second_output_file << acceptanceName << " " << observablename << " Unfolded: " << Afb << " +/-  " << AfbErr << endl;

        GetAfb(denominatorM, Afb, AfbErr);
        cout << " True Top from acceptance denominator: " << Afb << " +/-  " << AfbErr << "\n";
        second_output_file << acceptanceName << " " << observablename << " True_Top_from_acceptance_denominator: " << Afb << " +/-  " << AfbErr << "\n";

        Float_t Afb_true = Afb;

        TH1D *hPull = new TH1D ("pull", "pull", 100, -4, 4);


        for (int iPE = 0; iPE < lastPE; ++iPE)
        {
            GetCorrectedAfb(hPseudoData_unfolded[iPE], m_correctE_PEs[iPE], Afb, AfbErr);
            hPull->Fill( (Afb - Afb_true) / AfbErr );
        }

        TCanvas *c_pull = new TCanvas("c_pull", "c_pull", 500, 500);
        gStyle->SetOptStat(1111);
        hPull->Draw();
        c_pull->SaveAs("Pull_" + acceptanceName + Region + ".pdf");


        vector<double> afb_bins;
        vector<double> afb_bins;
        GetCorrectedAfbBinByBin(hData_unfolded, m_correctE, afb_bins, afb_bins, second_output_file);

        //scale to total xsec with option "width",  so that differential xsec is plotted
        hData_unfolded->Scale(1. / hData_unfolded->Integral(), "width");
        hTrue->Scale(1. / hTrue->Integral(), "width");


        for (int i = 1; i < nbins1D + 1; i++)
        {
            cout << i << " bin = " << hData_unfolded->GetBinContent(i) << " +/- " << hData_unfolded->GetBinError(i) << endl;
            second_output_file << acceptanceName << " " << observablename << " bin" << i << ": " << hData_unfolded->GetBinContent(i) << " +/- " << hData_unfolded->GetBinError(i) << endl;
        }


        TCanvas *c_test = new TCanvas("c_final", "c_final", 500, 500);

        hData_unfolded->GetXaxis()->SetTitle(xaxislabel);
        hData_unfolded->GetYaxis()->SetTitle("1/#sigma d#sigma/d(" + xaxislabel + ")");
        hData_unfolded->SetMinimum(0.0);
        hData_unfolded->SetMaximum( 2.0 * hData_unfolded->GetMaximum());
        hData_unfolded->SetMarkerStyle(23);
        hData_unfolded->SetMarkerSize(1.5);
        hData_unfolded->Draw("E");
        hData_unfolded->SetLineWidth(lineWidth);
        hTrue->SetLineWidth(lineWidth);
        hTrue->SetLineColor(TColor::GetColorDark(kGreen));
        hTrue->SetFillColor(TColor::GetColorDark(kGreen));
        hTrue->SetFillStyle(3353);
        hTrue->Draw("hist same");
        hData_unfolded->Draw("E same");

        TLegend *leg1 = new TLegend(0.55, 0.62, 0.9, 0.838, NULL, "brNDC");
        leg1->SetEntrySeparation(100);
        leg1->SetFillColor(0);
        leg1->SetLineColor(0);
        leg1->SetBorderSize(0);
        leg1->SetTextSize(0.03);

        leg1->AddEntry(hData_unfolded, "( Data - BG ) Unfolded");
        leg1->AddEntry(hTrue,    "SM parton level (mc@nlo)", "F");
        leg1->Draw();

        TPaveText *pt1 = new TPaveText(0.19, 0.85, 0.42, 0.89, "brNDC");
        pt1->SetName("pt1name");
        pt1->SetBorderSize(0);
        pt1->SetFillStyle(0);

        TText *blah;
        blah = pt1->AddText("CMS Preliminary, 5.0 fb^{-1} at  #sqrt{s}=7 TeV");
        blah->SetTextSize(0.032);
        blah->SetTextAlign(11);
        pt1->Draw();

        c_test->SaveAs("finalplot_unfolded_" + acceptanceName + Region + ".pdf");

        ch_data->Delete();

        ch_top->Delete();

        for (int iBkg = 0; iBkg < nBkg; ++iBkg)
        {
            ch_bkg[iBkg]->Delete();
        }

    }

    myfile.close();
    second_output_file.close();
}

#ifndef __CINT__
int main ()
{
    AfbUnfoldExample();    // Main program when run stand-alone
    return 0;
}
#endif

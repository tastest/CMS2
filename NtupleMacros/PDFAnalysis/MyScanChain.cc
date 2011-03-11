
//
// Dave "the one but not the only" Evans 
//

#include "MyScanChain.h"

// ROOT includes
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

// CMS2 includes
#include "../CORE/CMS2.h"
#include "../Tools/DileptonHypType.h"
#include "../Tools/tools.cc"
#include "../CORE/utilities.h"

// LHAPDF
#include "/tas07/disk00/jribnik/lhapdf/include/LHAPDF/LHAPDF.h"

//
// Namespaces
//
using namespace tas;

//
// names of analysis cuts
//

//
// for jets binning
//
enum JetBin {
    JET_0J,
    JET_1J,
    JET_GT2J,
    JET_ANY,
};

static const char jetbin_names[][128] = { "0j", "1j", "2j", "allj"};

//
// functions
//

bool MyScanChain::PassAnalysisSelection()
{

    //
    // does this event pass the analysis selection
    //

    // mc leptons
    std::vector<unsigned int> mcLeptonIndices;
    int nGoodLep = 0;
    for (size_t i = 0; i < cms2.genps_id().size(); ++i)
    {
        if (!(abs(cms2.genps_id()[i]) == 11 || abs(cms2.genps_id()[i]) == 13)) continue;
        if ( cms2.genps_p4()[i].Pt() < 20.0 || abs(cms2.genps_p4()[i].Eta()) > 2.5) continue;
        nGoodLep++;
        mcLeptonIndices.push_back(i);
    }

    // mc jets
    /*
    int nGoodJet = 0;
    for (size_t j = 0; j < cms2.evt_ngenjets(); ++j) 
    {
        if (cms2.genjets_p4()[j].Pt() < 30.0) continue;
        if (fabs(cms2.genjets_p4()[j].Eta()) > 2.5) continue;
        bool clean = true;
        for ( size_t i = 0; i < mcLeptonIndices.size(); ++i) 
        {
            if (ROOT::Math::VectorUtil::DeltaR(cms2.genjets_p4()[j], cms2.genps_p4()[mcLeptonIndices[i]]) < 0.4) {
                clean = false;
                break;
            }
        }
        if (clean) nGoodJet ++;
    }
    */

    //if (nGoodLep == 2 && nGoodJet >= 2) return = true;
    if (nGoodLep == 2) return true;
    return false;

}

float MyScanChain::GetGenMeff()
{
    //
    // compute MC quantities to plot
    //

    // get all the status 3 b quarks and leptons
    float scalarSum_g_pt = 0.0;
    float scalarSum_obslep_pt = 0.0;
    float scalarSum_quark_pt = 0.0; 
    LorentzVector vector_nu;

    unsigned int nLep = 0;

    for (size_t i = 0;  i < cms2.genps_id().size(); ++i) {

        if (cms2.genps_status()[i] == 3 && cms2.genps_id_mother()[i] != 21212) {
            int pdgId = abs(cms2.genps_id()[i]);
            //
            // with fiduciality cuts
            //
            // scalar sum of observable lepton pt
            if (pdgId == 11 || pdgId == 13) {
                if (fabs(cms2.genps_p4()[i].Eta()) < 2.5 && cms2.genps_p4()[i].Pt() > 20.0) {
                    scalarSum_obslep_pt += cms2.genps_p4()[i].Pt();
                    nLep ++;
                }
            }
            // quarks
            // don't count tops because they seem to decay at status 3
            else if (pdgId == 1 || pdgId == 2 || pdgId == 3 || pdgId == 4 || pdgId == 5) {
                if (fabs(cms2.genps_p4()[i].Eta()) < 2.5 && cms2.genps_p4()[i].Pt() > 30.0) {
                    scalarSum_quark_pt += cms2.genps_p4()[i].Pt();
                }
            }

            // gluons
            // assume that the gen particles will not contain
            // a gluon which is the daughter of a particle which is not a gluon
            // which is the daughter of a gluon.  Thus it is adequate to select
            // gluons whos parent is not a gluon
            else if (pdgId == 21) {
                if (fabs(cms2.genps_p4()[i].Eta()) < 2.5 && cms2.genps_p4()[i].Pt() > 20.0) {
                    scalarSum_g_pt += cms2.genps_p4()[i].Pt();
                }
            }

            // invisible
            else if (pdgId == 12 || pdgId == 14 || pdgId == 1000022) {
                // vector sum
                vector_nu += cms2.genps_p4()[i];
            }

        }

    }

    if (nLep > 1) return scalarSum_g_pt + scalarSum_obslep_pt + scalarSum_quark_pt + vector_nu.Pt();
    else return 0.0;

}

void MyScanChain::Fill(TH1F** hist, const unsigned int hyp, const float &val, const float &weight)
{
    hist[hyp]->Fill(val, weight);
    hist[DILEPTON_ALL]->Fill(val, weight);
}

void MyScanChain::FormatHist(TH1F** hist, std::string sampleName, std::string name, int n, float min, float max)
{       
    // loop on EB, EE
    for (unsigned int i = 0; i < 4; ++i)
    {
        std::string str = dilepton_hypo_names[i];
        std::string title = name + "_" + str;
        hist[i] = new TH1F(Form("%s_%s_%s", sampleName.c_str(), name.c_str(), str.c_str()),
                title.c_str(), n, min, max);
        hist[i]->GetXaxis()->SetTitle(name.c_str());
        hist[i]->Sumw2();
    }
}    

//
// set the pdf to use
//

void MyScanChain::specifyPDF(std::string pdfName1, std::string pdfName2) 
{
    pdfName1_ = pdfName1;
    pdfName2_ = pdfName2;
}

//
// Main function
//

int MyScanChain::ScanChain(std::string sampleName, TChain *chain, float kFactor, int nEvents, std::string skimFilePrefix) {

    std::cout << "scanning " << sampleName << std::endl;
    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    if (rootdir == 0){
        std::cout<<"Head directory root: not found. Try Rint: ..."<<std::endl;
        rootdir = gROOT->GetDirectory("Rint:");
        if (rootdir){
            std::cout<<"OK: Got Rint:"<<std::endl;
        } else {
            std::cout<<"ERROR: no root: or Rint: found. Histograms will likely be lost"<<std::endl;
        }
    } 
    rootdir->cd(); 

    //
    // set up LHAPDF
    //

    // initialise the pdf 
    const int pdfSubset = 0;
    LHAPDF::initPDFSet(1, pdfName1_, LHAPDF::LHGRID, pdfSubset);
    LHAPDF::initPDFSet(2, pdfName2_, LHAPDF::LHGRID, pdfSubset);

    // the weight for the central value for 
    // some other pdf
    double pdfWeightOther = 1.0;

    // set up array to hold the weights for each event
    // for the pdf used to generate the sample
    // size is set to maximum possible
    double pdfWeights[MAXWEIGHT];
    unsigned int nsets = 1;
    int genPdf = 1;
    int otherPdf = 2;
    if (LHAPDF::numberPDF(genPdf) > 1) nsets += LHAPDF::numberPDF(genPdf);

    //
    // format histograms
    // and set up quantities to study
    // for the pdf used to generate sample
    //

    Int_t nbins_meff = 2000;
    Float_t min_meff = 0.0;
    Float_t max_meff = 2000.0;
    TH1F   *histArr[MAXWEIGHT];
    double acceptance[MAXWEIGHT];
    double nTotal[MAXWEIGHT];
    double nPass[MAXWEIGHT];
    for (unsigned int i = 0; i < MAXWEIGHT; ++i) {
        acceptance[i] = 0.0;
        nTotal[i] = 0.0;
        nPass[i] = 0.0;
        histArr[i] = new TH1F(Form("h1_pass_%i", i), Form("pass_%i", i), nbins_meff, min_meff, max_meff);
    }

    TH1F *h1_centre   = new TH1F(Form("%s_h1_centre_all", sampleName.c_str()),  "centre", nbins_meff, min_meff, max_meff);
    TH1F *h1_up   = new TH1F(Form("%s_h1_up_all", sampleName.c_str()),  "up", nbins_meff, min_meff, max_meff);
    TH1F *h1_down  = new TH1F(Form("%s_h1_down_all", sampleName.c_str()), "down", nbins_meff, min_meff, max_meff);

    //
    // and for the other specified pdf
    //
    double acceptanceOther = 0.0;
    double nTotalOther = 0.0;
    double nPassOther = 0.0;

    //
    // file loop
    //

    unsigned int nEventsChain=0;
    if(nEvents == -1) nEvents = chain->GetEntries();
    nEventsChain = nEvents;
    unsigned int nEventsTotal = 0;
    int i_permille_old = 0;

    TObjArray *listOfFiles = chain->GetListOfFiles();
    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;
    while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {
        TFile *f = TFile::Open(currentFile->GetTitle());
        TTree *tree = (TTree*)f->Get("Events");
        cms2.Init(tree);

        //Event Loop
        ULong64_t nEvents = tree->GetEntries();
        for(ULong64_t event = 0; event < nEvents; ++event) {
            cms2.GetEntry(event);
            ++nEventsTotal;

            // Progress feedback to the user
            int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
            if (i_permille != i_permille_old) {
                // xterm magic from L. Vacavant and A. Cerri
                if (isatty(1)) {
                    printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                            "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
                    fflush(stdout);
                }
                i_permille_old = i_permille;
            }

            // work out event weight
            float weight = cms2.evt_scale1fb() * kFactor;

            //
            // deal with DY
            //

            if ((sampleName == "dyee" || sampleName == "dymm" || sampleName == "dytt")) {

                // if we are not looking at the madgraph then we must have one of the pythia pieces
                // - if these have an event with mass above the start of the madgraph sample then skip

                if (TString(cms2.evt_dataset()).Contains("madgraph") == false) {
                    bool skipEvent = false;
                    for (unsigned int i = 0; i < genps_p4().size(); i++){
                        if(abs(cms2.genps_id()[i]) == 23 && cms2.genps_p4()[i].M() > 50.) {
                            skipEvent = true;
                            break;
                        }
                    }

                    // skip the event if necessary
                    if (skipEvent) continue;

                    // if the event was not skipped then set the k factor
                    // for the appropriate sample

                    if (TString(evt_dataset()).Contains("M10to20") == true) { //10 < mll < 20 
                        weight = weight*3457./2659.;
                    } else { // 20 < mll < 50
                        weight = weight * 1666./1300.;
                    }

                    // otherwise we are looking at the madgraph sample
                    // this has a mass > 50.0 and we always include it

                } else {
                    // set the weight for madgraph
                    weight = weight*3048./2400.;
                }

            }

            //
            // do PDF analysis
            //

            // central value for the weight for the pdf used to generate the sample
            // is always going to be one
            pdfWeights[0] = 1.0;

            // calculate the central value of the other pdf specified for this event
            LHAPDF::usePDFMember(otherPdf, 0);
            double fx1Q0Other = LHAPDF::xfx(cms2.pdfinfo_x1(), cms2.pdfinfo_scale(), cms2.pdfinfo_id1()) / cms2.pdfinfo_x1();
            double fx2Q0Other = LHAPDF::xfx(cms2.pdfinfo_x2(), cms2.pdfinfo_scale(), cms2.pdfinfo_id2()) / cms2.pdfinfo_x2();

            // calculate the central value of the pdf used to generate the sample
            LHAPDF::usePDFMember(genPdf, 0);
            double fx1Q0 = LHAPDF::xfx(cms2.pdfinfo_x1(), cms2.pdfinfo_scale(), cms2.pdfinfo_id1()) / cms2.pdfinfo_x1();
            double fx2Q0 = LHAPDF::xfx(cms2.pdfinfo_x2(), cms2.pdfinfo_scale(), cms2.pdfinfo_id2()) / cms2.pdfinfo_x2();

            // calculate a weight for this events central value into the other pdf
            pdfWeightOther = (fx1Q0Other*fx2Q0Other)/(fx1Q0*fx2Q0);

            // calculate the weight for the ith subset
            // for this event
            for (unsigned int subset = 1; subset < nsets; ++subset) 
            {
                LHAPDF::initPDF(subset);
                double fx1Qi = LHAPDF::xfx(cms2.pdfinfo_x1(), cms2.pdfinfo_scale(), cms2.pdfinfo_id1()) / cms2.pdfinfo_x1();
                double fx2Qi = LHAPDF::xfx(cms2.pdfinfo_x2(), cms2.pdfinfo_scale(), cms2.pdfinfo_id2()) / cms2.pdfinfo_x2();
                pdfWeights[subset] = (fx1Qi*fx2Qi)/(fx1Q0*fx2Q0);
            }

            //
            // histogram the observable quantities for each eigenvector variation weight
            // and record the acceptance for passing the analysis selection for each eigenvector variation weight
            //

            bool passAnalysisSelection = PassAnalysisSelection();
            float genmeff = GetGenMeff();

            for (unsigned int subset = 0; subset < nsets; ++subset) 
            {

                // acceptance calculation
                nTotal[subset] += pdfWeights[subset];
                nTotalOther += pdfWeightOther;
                if (passAnalysisSelection) {
                    nPass[subset] += pdfWeights[subset];
                    nPassOther += pdfWeightOther;
                }

                // distribution calculations
                // in this case it's the meff w.r.t. events passing the analysis selection
                if (passAnalysisSelection) {
                    histArr[subset]->Fill(genmeff, pdfWeights[subset]*weight);
                }

            }

        } // end loop on files

    } // end loop on events

    if ( nEventsChain != nEventsTotal ) {
        std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
    }

    //
    // Analyse the observables
    //

    //
    // calculate the weighted acceptance for each eigenvector variation
    //

    acceptanceOther = nPassOther / nTotalOther;
    for (unsigned int i = 0; i < nsets; ++i) {
        acceptance[i] = nPass[i] / nTotal[i];
    } 

    //
    // calculate the variation in the acceptance from the maximum variation
    // of all the eigenvector variations
    //

    double plus_max = 0.0;
    double minus_max = 0.0;
    double X0 = acceptance[0];

    for (unsigned int subset = 0; subset < ((nsets - 1)/2); ++subset)
    {
        double Xi_up = acceptance[(subset*2) + 1];
        double Xi_down = acceptance[(subset*2) + 2];
        plus_max += pow(max(max(Xi_up - X0, Xi_down - X0), 0.0), 2);
        minus_max += pow(max(max(X0 - Xi_up, X0 - Xi_down), 0.0), 2);
    }

    plus_max = sqrt(plus_max);
    minus_max = sqrt(minus_max);

    //
    // now do the same for the quantities as a function of some variable
    // (histograms)
    //

    for (Int_t bin = 1; bin <= nbins_meff; ++bin)
    {

        Double_t plus_max = 0;
        Double_t minus_max = 0;
        Double_t X0 = histArr[0]->GetBinContent(bin);
        Double_t binCenter = histArr[0]->GetBinCenter(bin);

        if (X0 == 0) continue;

        for (unsigned int subset = 0; subset < ((nsets - 1)/2); ++subset)
        {
            Double_t Xi_up = histArr[(subset*2) + 1]->GetBinContent(bin);
            Double_t Xi_down = histArr[(subset*2) + 2]->GetBinContent(bin);
            plus_max += pow(max(max(Xi_up - X0, Xi_down - X0), 0.0), 2);
            minus_max += pow(max(max(X0 - Xi_up, X0 - Xi_down), 0.0), 2);
        }

        plus_max = sqrt(plus_max);
        minus_max = sqrt(minus_max);

        h1_centre->SetBinContent(bin, X0);
        h1_down->SetBinContent(bin, X0 - minus_max);
        h1_up->SetBinContent(bin, X0 + plus_max);

    }

    //
    // print out the acceptance results
    // (the results for the histograms are written to a root file)
    //

    std::cout << "[MyScanChain::ScanChain] " << sampleName << std::endl;
    std::cout << "[MyScanChain::ScanChain] Analysing PDF uncertainty on the acceptance" << std::endl;
    std::cout << "[MyScanChain::ScanChain] Central value is         : " << X0 << "$^{+" << plus_max << "}_{-" << minus_max << "}" << std::endl;
    std::cout << "[MyScanChain::ScanChain] Relative uncertainty is  : $^{+" << plus_max/X0 << "}_{-" << minus_max/X0 << "}" << std::endl;
    std::cout << "[MyScanChain::ScanChain] Central value for " << pdfName2_ << ": " << acceptanceOther << std::endl;

    //
    // clean up the temporary hists
    //

    for (unsigned int i = 0; i < MAXWEIGHT; ++i) delete histArr[i];

    //
    // make sure we're back in the right root dir
    //

    rootdir = gROOT->GetDirectory("root:");
    if (rootdir) rootdir->cd();
    else{
        std::cout<<"Cant find root: . Current dir is "<<gDirectory->GetName()<<std::endl;
        rootdir = gROOT->GetDirectory("Rint:");
        if (rootdir){
            std::cout<<"OK, got Rint: "<<std::endl;
            rootdir->cd();
        } else {
            std::cout<<"Cant find Rint: either . Current dir is "<<gDirectory->GetName()<<std::endl;
        }
    }

    return 0;
}


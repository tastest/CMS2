
//
// ttbar -> ll
// Dave "the one but not the only" Evans 
//

#include "MyScanChain.h"

// ROOT includes
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"

#include "Math/LorentzVector.h"

// CMS2 includes
#include "CMS2.h"
#include "../../CORE/trackSelections.h"
#include "../../CORE/electronSelections.h"
#include "../../Tools/DileptonHypType.h"

//
// Namespaces
//
using namespace tas;

//
//
//

enum ele_selection {
    PASS_PT,
    PASS_NOSECOND,
    PASS_ISFIDUCIAL,
    PASS_MET,
    PASS_ISO,
    PASS_JETVETO,

    PASS_DPHI,
    PASS_DETA,
    PASS_HOE,
    PASS_LSHAPE,
    PASS_D0,
    PASS_DETA_CAND02,
    PASS_LSHAPE_CAND02,
    PASS_EXTRA,
    PASS_CONV,
    PASS_NOMUON,

    PASS_CLEANEVENT,

};

double dRbetweenVectors(const LorentzVector &vec1, const LorentzVector &vec2 ){

    double dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
    double deta = vec1.Eta() - vec2.Eta();
    return sqrt(dphi*dphi + deta*deta);
} 

// to decdie if to fill the EB histo (zero in the array)
// or the EE histo (one in the array)
enum DetType MyScanChain::getSubdet(int eleIndex)
{
    // determine what detector the electron is in
    if (cms2.els_fiduciality()[eleIndex] & (1<<ISEB)) return DET_EB;
    else if (cms2.els_fiduciality()[eleIndex] & (1<<ISEE)) return DET_EE;
    std::cout << "ERROR! not in EB or EE" << std::endl;
    return DET_ALL;
}

void MyScanChain::Format2DHist(TH2F** hist, std::string name, Int_t nx, Float_t minx, Float_t maxx,
        Int_t ny, Float_t miny, Float_t maxy)
{
    // loop on EB, EE
    for (unsigned int i = 0; i < 3; ++i)
    {
        std::string det = det_names[i];
        hist[i] = new TH2F(Form("%s_%s_%s", sampleName_.c_str(), name.c_str(), det.c_str()),
                name.c_str(), nx, minx, maxx, ny, miny, maxy);
    }
}

void MyScanChain::FormatHist(TH1F** hist, std::string name, Int_t n, Float_t min, Float_t max)
{
    // loop on EB, EE
    for (unsigned int i = 0; i < 3; ++i)
    {
        std::string det = det_names[i];
        hist[i] = new TH1F(Form("%s_%s_%s", sampleName_.c_str(), name.c_str(), det.c_str()),
                name.c_str(), n, min, max);
        //hist[i]->SetFillColor(sampleName_.histo_color);
        hist[i]->GetXaxis()->SetTitle(name.c_str());
    }
}

void MyScanChain::FormatEffHist(EffMulti** hist, bool lessThan, float thresholdEB, float thresholdEE, std::string name)
{
    // loop on EB, EE
    for (unsigned int i = 0; i < 3; ++i)
    {
        std::string det = det_names[i];
        hist[i] = new EffMulti(lessThan, thresholdEB, thresholdEE, sampleName_, name, det);
    }
}

bool MyScanChain::CheckCutsNM1(cuts_t apply, cuts_t remove, cuts_t passed)
{           
    if ((passed & (apply & (~remove))) == (apply & (~remove))) return true;
    return false;
}   

bool MyScanChain::CheckCuts(cuts_t apply, cuts_t passed)
{
    if ((apply & passed) == apply) return true;
    return false;
}


//
// Main function
//
int MyScanChain::ScanChain(bool isData, std::string sampleName, TChain *chain, int nEvents, std::string skimFilePrefix) {

    //
    // define counters
    //

    // count the (weighted and unweighted) number of candidates passing our cuts
    float             cands_passing[3];
    float             cands_passing_w2[3];
    unsigned int       cands_count[3];
    memset(cands_passing   , 0, sizeof(cands_passing       ));
    memset(cands_passing_w2        , 0, sizeof(cands_passing_w2    ));
    memset(cands_count             , 0, sizeof(cands_count         ));

    // set sampleName
    sampleName_ = sampleName;

    //
    //
    //
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

    //
    // format histograms
    //

    // General
    //
    FormatHist(h1_pt_, "pt", 100, 0, 100);
    FormatHist(h1_eta_, "eta", 100, -3, 3);
    FormatHist(h1_phi_, "phi", 100, -4, 4);

    // N-1
    FormatHist(h1_nm1_tcmet_, "nm1_tcmet", 20, 0, 100);
    FormatHist(h1_nm1_pfmet_, "nm1_pfmet", 20, 0, 100);
    FormatHist(h1_nm1_jetveto_, "nm1_jetveto", 100, 0, 100);
    FormatHist(h1_nm1_iso_, "nm1_iso", 100, 0, 1);
    FormatHist(h1_nm1_secondpt_, "nm1_secondpt", 100, 0, 100);

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

            //			float weight = cms2.evt_scale1fb()*0.01;
            float weight = 1.0;

            //
            // do event cleaning
            //
            // require bit 40 or 41 passed
            //
            if (!(cms2.l1_techbits2() & (1<<8) || cms2.l1_techbits2() & (1<<9))) continue;
            // require the anti-coincidence of 42 and 43 DIDN't pass
            //
            if (!(cms2.l1_techbits2() & (1<<10)) && (cms2.l1_techbits2() & (1<<11))) continue;
            if ((cms2.l1_techbits2() & (1<<10)) && !(cms2.l1_techbits2() & (1<<11))) continue;
            // require bits 36-39 DIDN't pass
            //
            if (cms2.l1_techbits2() & (1<<7) || cms2.l1_techbits2() & (1<<6) || 
                    cms2.l1_techbits2() & (1<<5) || cms2.l1_techbits2() & (1<<4)) continue;
            //
            // require bit zero for beams 
            if (isData && !(cms2.l1_techbits1() & (1<<0))) continue;
            //
            // at least 1 good vertex
            // count good vertexs
            //              
            int nGoodVertex = 0;
            for (size_t v = 0; v < cms2.vtxs_position().size(); ++v) {

                if (cms2.vtxs_isFake()[v]) continue;
                if (cms2.vtxs_ndof()[v] < 4) continue;
                if (fabs(cms2.vtxs_position()[v].Z()) > 15.0) continue;
                nGoodVertex ++;
            }
            if (nGoodVertex == 0) continue;

            //
            // if >= 10 tracks, require at least 25% high purity
            if (cms2.trks_ndof().size() >= 10) {
                int nHighPurityTracks = 0;
                for (size_t t = 0; t < cms2.trks_ndof().size(); ++t)
                    if (isTrackQuality(t, (1<<highPurity))) nHighPurityTracks ++;
                if (cms2.trks_ndof().size()/float(nHighPurityTracks) < 0.25) continue;
            }

            //
            // end event cleaning
            //



            // find candidate electron
            // assumes electrons are sorted by pT descending
            int eleIndex = 0;
            int eleSecondIndex = 0;
            bool foundFirst = false;
            bool foundSecond = false;
            for (size_t i = 0; i < cms2.evt_nels(); ++i)
            {       
                // no particle flow
                if (! cms2.els_type()[i] & (1<<ISECALDRIVEN)) continue;

                if (foundFirst && !foundSecond) {
                    eleSecondIndex = i;
                    foundSecond = true;
                    break;
                }       
                if (!foundFirst) { 
                    eleIndex = i;
                    foundFirst = true;
                }
            }

            // must have found first electron
            if (!foundFirst) continue;

            // get subdetector (for histogramming)
            int det = getSubdet(eleIndex);

            //
            // work out what cuts this event passes
            //

            // the cuts that this event passes
            cuts_t cuts_passed = 0;

            // pt cut
            if (cms2.els_p4()[eleIndex].Pt() > 20.0) cuts_passed |= (1<<PASS_PT);

            // don't allow events with a second electron above 20.0 GeV
            float secondPt = 0.0;
            if (foundSecond) secondPt == cms2.els_p4()[eleSecondIndex].Pt();
            if (secondPt < 20.0) cuts_passed |= (1<<PASS_NOSECOND);

            // impose fiducial cuts in Eta
            if (fabs(cms2.els_etaSC()[eleIndex]) < 1.4442 
                    || (fabs(cms2.els_etaSC()[eleIndex]) > 1.560 && fabs(cms2.els_etaSC()[eleIndex]) < 2.500)) cuts_passed |= (1<<PASS_ISFIDUCIAL);

            // isolation
            //std::cout << cms2.els_tkJuraIso().at(eleIndex) << std::endl;
            //			float iso_relsusy = electronIsolation_relsusy_cand1(eleIndex, true);


            //printf("cms2 = 0x%x\n", &cms2);
            float sum = cms2.els_tkJuraIso().at(eleIndex);
            if (fabs(cms2.els_etaSC().at(eleIndex)) > 1.479) sum += cms2.els_ecalIso().at(eleIndex);
            if (fabs(cms2.els_etaSC().at(eleIndex)) <= 1.479) sum += max(0., (cms2.els_ecalIso().at(eleIndex) -1.));
            sum += cms2.els_hcalIso().at(eleIndex);
            double pt = cms2.els_p4().at(eleIndex).pt();
            float iso_relsusy =  sum/max(pt, 20.);

            if (iso_relsusy < 0.10) cuts_passed |= (1<<PASS_ISO);

            // met
            if (cms2.evt_tcmet() > 30.0) cuts_passed |= (1<<PASS_MET);

            // jet veto... leading pT JPT jet that is dR > 0.4 from the nearest electron
            float leadingJPT = 0.0;
            int leadingJPTIndex = 0;
            for (size_t j = 0; j < cms2.jpts_p4().size(); ++j)
            {
                if ( TMath::Abs(cms2.jpts_p4()[j].eta()) > 2.5 ) continue;
                if ( TMath::Abs(dRbetweenVectors(cms2.els_p4()[eleIndex], cms2.jpts_p4()[j])) < 0.4) continue;

                // find leading pT JPT
                if (cms2.jpts_p4()[j].Et() > leadingJPT) {
                    leadingJPT = cms2.jpts_p4()[j].Et();
                    leadingJPTIndex = j;
                }
            }
            if (leadingJPT < 30.0) cuts_passed |= (1<<PASS_JETVETO);

            //
            // do plotting
            //


            h1_pt_[det]->Fill(cms2.els_p4()[eleIndex].Pt(), weight);
            h1_pt_[DET_ALL]->Fill(cms2.els_p4()[eleIndex].Pt(), weight);			

            h1_eta_[det]->Fill(cms2.els_etaSC()[eleIndex], weight);
            h1_eta_[DET_ALL]->Fill(cms2.els_etaSC()[eleIndex], weight);

            h1_phi_[det]->Fill(cms2.els_phiSC()[eleIndex], weight);
            h1_phi_[DET_ALL]->Fill(cms2.els_phiSC()[eleIndex], weight);

            const cuts_t pass_all = (1<<PASS_PT) | (1<<PASS_NOSECOND) | (1<<PASS_ISFIDUCIAL) | (1<<PASS_ISO) | (1<<PASS_MET) | (1<<PASS_JETVETO);

            if (CheckCutsNM1(pass_all, (1<<PASS_MET), cuts_passed)) {
                h1_nm1_tcmet_[det]->Fill(cms2.evt_tcmet(), weight);
                h1_nm1_tcmet_[DET_ALL]->Fill(cms2.evt_tcmet(), weight);

                h1_nm1_pfmet_[det]->Fill(cms2.evt_pfmet(), weight);
                h1_nm1_pfmet_[DET_ALL]->Fill(cms2.evt_pfmet(), weight);

            }

            if (CheckCutsNM1(pass_all, (1<<PASS_JETVETO), cuts_passed)) {
                h1_nm1_jetveto_[det]->Fill(leadingJPT, weight);
                h1_nm1_jetveto_[DET_ALL]->Fill(leadingJPT, weight);
            }

            if (CheckCutsNM1(pass_all, (1<<PASS_ISO), cuts_passed)) {
                h1_nm1_iso_[det]->Fill(iso_relsusy, weight);
                h1_nm1_iso_[DET_ALL]->Fill(iso_relsusy, weight);
            }

            if (CheckCutsNM1(pass_all, (1<<PASS_NOSECOND), cuts_passed)) {
                h1_nm1_secondpt_[det]->Fill(secondPt, weight);
                h1_nm1_secondpt_[DET_ALL]->Fill(secondPt, weight);
            }

            //
            // Count
            //

            cands_passing[det] += weight;
            cands_passing_w2[det] += weight * weight;
            cands_count[det] ++;

            cands_passing[DET_ALL] += weight;
            cands_passing_w2[DET_ALL] += weight * weight;
            cands_count[DET_ALL] ++;

        } // end loop on files

    } // end loop on events

    if ( nEventsChain != nEventsTotal ) {
        std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
    }

    //
    // print table entry
    //

    std::cout.flush();
    std::cout << std::endl;
    for (unsigned int i = 0; i < 3; ++i) {
        std::string str = det_names[i];
        std::cout << " & " << det_names[i] << "\t";
    }
    std::cout << "\\\\ \\hline" << std::endl;
    std::cout << sampleName << "\t";
    for (unsigned int i = 0; i < 3; ++i) {
        std::cout << " & " << cands_passing[i] << " $\\pm$ " << sqrt(cands_passing_w2[i]) << "\t";
    }
    std::cout << "\\\\ \\hline" << std::endl;

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


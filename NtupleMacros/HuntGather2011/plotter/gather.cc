#include "BabySample.h"
#include "cuts.h"
#include "gather.h"
#include "../../Tools/goodrun.h"
#include "BabyDorkIdentifier.h"

#include "TCanvas.h"
#include "TChain.h"
#include "TCut.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPRegexp.h"
#include "TROOT.h"
#include "TTreePlayer.h"

#include <algorithm>
#include <vector>
#include <iostream>

float GetIntLumi(TChain *c, float lumi, int brun, int bls, int erun, int els)
{

    TCut c_goodrun(Form("((run>%i&&run<%i)||(run==%i&&ls>=%i)||(run==%i&&ls<=%i))&&goodrun_json(run,ls)", brun, erun, brun, bls, erun, els));

    // goodrun plus events beyond range of goodrun
    // which are not goodrun penalized
    TCut c_goodrunplus(Form("(((run>%i&&run<%i)||(run==%i&&ls>=%i)||(run==%i&&ls<=%i))&&goodrun_json(run,ls))||(run>%i||(run==%i&&ls>%i))", brun, erun, brun, bls, erun, els, erun, erun, els));

    // brun:bls -> erun:els
    reset_babydorkidentifier();
    int n_goodrun = c->GetEntries(c_goodrun+inclusivez_dilep); 
    // total
    reset_babydorkidentifier();
    int n_total   = c->GetEntries(c_goodrunplus+inclusivez_dilep);
    // that which is new
    int n_new     = n_total-n_goodrun;

    float newlumi = ((float)(n_new*lumi))/(float)n_goodrun;

    std::cout << "in the good run list for " << lumi << " there are " << n_goodrun << std::endl;
    std::cout << "out of the good run list there are " << n_new << std::endl;
    std::cout << "that means a new lumi of " << newlumi << std::endl;

    return lumi+newlumi;
}

float GetIntLumi(TChain *c, float lumi)
{
    int brun = min_run();
    int bls  = min_run_min_lumi();
    int erun = max_run();
    int els  = max_run_max_lumi();
    std::cout << brun << ", " << bls << " : " << erun << ", " << els << std::endl;
    return GetIntLumi(c, lumi, brun, bls, erun, els);
}

TCanvas* DrawAll(TCut var, const char *savename, TCut sel, TCut presel, float intlumipb, unsigned int nbins, float xlo, float xhi, bool integrated,
        std::vector<BabySample*> bss)
{

    std::vector<TH1F*> vh_background;
    std::vector<TH1F*> vh_signal;
    std::vector<TH1F*> vh_data;
    TH1F* buffer;

    //
    // loop on baby samples
    //

    for(unsigned int i = 0; i < bss.size(); ++i) {

        // construct the selection for this sample
        // and get the plot of that selection
        TCut selection = sel + presel;
        buffer = Plot (bss[i], var, selection, intlumipb, nbins, xlo, xhi, integrated, gDrawAllCount);

        // if the plot doesn't already exist then
        // add it to the appropriate vector of plots
        if (bss[i]->type() == DATA && (find(vh_data.begin(), vh_data.end(), buffer) == vh_data.end())) {
            vh_data.push_back(buffer);
        } else if (bss[i]->type() == BACKGROUND && find(vh_background.begin(), vh_background.end(), buffer) == vh_background.end()) {
            vh_background.push_back(buffer);
        } else if (bss[i]->type() == SIGNAL && find(vh_signal.begin(), vh_signal.end(), buffer) == vh_signal.end()) {
            vh_signal.push_back(buffer);
        }


    }

    // MUST have at least one data and at least one MC
    if (vh_background.size() == 0 || vh_data.size() == 0) {
        std::cout << "[DrawAll] ERROR - MUST have at least one MC and at least one data" << std::endl;
        return 0;
    }

    //
    // sort histograms by their total contributions
    // makes for prettier stack plots
    //

    sort(vh_background.begin(), vh_background.end(), sortHistsByIntegral);
    sort(vh_signal.begin(), vh_signal.end(), sortHistsByIntegral);
    sort(vh_data.begin(), vh_data.end(), sortHistsByIntegral);

    //
    // do the stacking
    // NOTE - signals are overlaid NOT stacked
    //

    makeStack(vh_background);
    makeStack(vh_data);

    //
    // set up the legend
    //

    TLegend* leg = new TLegend(0.7,0.5,0.95,0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetShadowColor(0);

    //
    // do the drawing
    //

    TCanvas *c1 = new TCanvas(savename);
    c1->SetTopMargin(0.08);
    c1->cd();

    // do the background MC histograms
    for(unsigned int i = 0; i < vh_background.size(); ++i) 
    { 
        // the zeroth sets the axis
        if (i == 0) {
            vh_background[0]->Draw("hist");
            vh_background[0]->GetXaxis()->SetNdivisions(504);
        }
        // the others are drawn the same
        else vh_background[i]->Draw("histsame");
        // add each histogram to the legend
        addToLegend(leg, vh_background.at(i), "f");
    }

    // do the signal MC histograms
    for(unsigned int i = 0; i < vh_signal.size(); ++i)
    {
        vh_signal[i]->Draw("histsame");
        addToLegend(leg, vh_signal.at(i), "l");
    }

    // do the data histogram
    vh_data[0]->Draw("samee1");
    addToLegend(leg, vh_data[0], "lp");

    // set the optimum y axis scale
    float ymax = vh_data[0]->GetMaximum() > vh_background[0]->GetMaximum() 
        ? vh_data[0]->GetMaximum()+2*sqrt(vh_data[0]->GetMaximum()) : vh_background[0]->GetMaximum()*1.25;
    vh_background[0]->SetMaximum(ymax);

    // draw the legend and tidy up
    leg->Draw();
    c1->RedrawAxis();
    gDrawAllCount++;

    return c1;
}

TH1F* Plot(const BabySample *bs, TCut var, TCut selection, float intlumipb, 
        unsigned int nbins, float xlo, float xhi, bool integrated, unsigned int isfx)
{

    //
    // Construct the name and title
    // and command to draw
    // 
    
    char *name = 0;
    if (integrated) name = Form("%s_%s_%s_%i_int", bs->pfx(), selection.GetName(), var.GetName(), gDrawAllCount);
    else name = Form("%s_%s_%s_%i", bs->pfx(), selection.GetName(), var.GetName(), gDrawAllCount);
    char *title = Form("%s_%s_%s, ~%.2f/pb", bs->pfx(), selection.GetName(), var.GetName(), 1e3*intlumipb);
    char *drawcommand = Form("%s>>+%s", var.GetTitle(), name);

    // set up the good run list
    int brun = min_run();
    int bls  = min_run_min_lumi();
    int erun = max_run();
    int els  = max_run_max_lumi();
    TCut c_goodrunplus(Form("(((run>%i&&run<%i)||(run==%i&&ls>=%i)||(run==%i&&ls<=%i))&&goodrun(run,ls))||(run>%i||(run==%i&&ls>%i))", 
                brun, erun, brun, bls, erun, els, erun, erun, els));

    if (!strcmp("CUT", var.GetName()))
        var.SetName(var.GetTitle());

    // construct the duplicate removal cut according
    // to the type of baby hypothesis
    TCut c_notduplicate = "";
    if (bs->chain()->GetBranch("pt3"))        c_notduplicate = ("! is_duplicate(run,evt,ls,pt1,pt2,pt3)");
    else if (bs->chain()->GetBranch("pt2"))   c_notduplicate = ("! is_duplicate(run,evt,ls,pt1,pt2)");
    else                                c_notduplicate = ("! is_duplicate(run,evt,ls,pt1)");

    // assemble the cut
    if (bs->type() == DATA) {
        selection += c_goodrunplus + c_notduplicate;
        std::cout << selection << std::endl;
    }

    //
    // set the normalisation scale
    //
    TCut scale = Form("scale1fb*(%f/1000.0)*%f", intlumipb, bs->kfactor());
    if (bs->type() == DATA) scale = "1.0";

    //
    // If the histogram does not already exist
    // then create it
    //

    TH1F *h = 0;
    if (!(h = (TH1F*)gROOT->FindObjectAny(name))) h = new TH1F(name, title, nbins, xlo, xhi);

    //
    // Draw the histogram
    //

    // fill an event list for the plot selection
    // and then get this event list
    bs->chain()->Draw(">>evtlist", selection, "goff");
    TEventList* elist = (TEventList*)gDirectory->Get("evtlist");

    // set the event list for events passing
    // the plot selection and draw the plot
    bs->chain()->SetEventList(elist);
    bs->chain()->Draw(drawcommand, scale, "goff");

    // if the events are data then
    // use the event list to print a scan
    if (bs->type() == DATA) {
        TTreePlayer *tp = (TTreePlayer*)bs->chain()->GetPlayer();
        tp->SetScanRedirect(kTRUE);
        tp->SetScanFileName(Form("%s_%s_%i_%s.out", selection.GetName(), var.GetName(), bs->pfx(), bs->pfx2()));
        // trileptons
        if (bs->chain()->GetBranch("pt3"))
            bs->chain()->Scan("dataset:run:ls:evt:pt1:pt2:pt3", "", "colsize=20");
        // dileptons
        else if (bs->chain()->GetBranch("pt2")) // dilep
            bs->chain()->Scan("dataset:run:ls:evt:pt1:pt2", "", "colsize=20");
        // single leptons
        else
            bs->chain()->Scan("dataset:run:ls:evt:pt1", "", "colsize=20");
    }

    // reset to an empty list
    TEventList* emptylist = 0;
    bs->chain()->SetEventList(emptylist);

    //
    // Set style options
    //

    if (bs->type() == DATA) {
        h->SetMarkerColor(bs->color());
        h->SetMarkerStyle(bs->style());
    } else if (bs->type() == BACKGROUND) {
        h->SetFillColor(bs->color());
        h->SetFillStyle(bs->style());
    } else if (bs->type() == SIGNAL) {
        h->SetLineColor(bs->color());
        h->SetLineStyle(bs->style());
        h->SetLineWidth(2);
    }

    h->SetTitle(Form("%s, ~%.2f/pb", selection.GetName(), intlumipb));
    if (integrated) {
        h = slideIntegrated(h);
        h->GetXaxis()->SetTitle(Form("integrated %s", var.GetName()));
    }
    else
        h->GetXaxis()->SetTitle(var.GetName());

    //
    // Move overflow to the last bin
    //

    float overflow = h->GetBinContent(nbins+1);
    h->SetBinContent(nbins,overflow);
    h->SetBinContent(nbins+1,0.);

    return h;

}

bool sortHistsByIntegral(TH1* h1, TH1* h2)
{
    return h1->Integral() > h2->Integral();
}

void addToLegend(TLegend *leg, TH1F *hist, TString opt)
{

    // For extracting sample prefix
    TPRegexp preg("^([^_]+)_.*$");
    TString  s_pfx("");

    s_pfx   = ((TObjString*)(preg.MatchS(TString(hist->GetName()))->At(1)))->GetString();
    leg->AddEntry(hist, s_pfx.Data(), opt);

}

void makeStack(std::vector<TH1F*> &v_hists)
{

    for(unsigned int i = 0; i < v_hists.size(); ++i) {
        for(unsigned int j = i+1; j < v_hists.size(); ++j) {
            (v_hists[i])->Add(v_hists[j], 1);
        }
    }

}


TH1F* slideIntegrated(TH1F* h1)
{
    int NbinsX = h1->GetNbinsX();
    for(int i = 1; i <= NbinsX; ++i) {
        // don't forget the overflow!
        float integral = h1->Integral(i,NbinsX+1);
        h1->SetBinContent(i,integral);
        h1->SetBinError(i,TMath::Sqrt(integral));
    }

    return h1;
}


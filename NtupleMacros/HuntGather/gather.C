#include "BabyDorkIdentifier.h"
#include "BabySample.h"
#include "cuts.h"
#include "gather.h"
#include "goodrun.h"

#include "TCanvas.h"
#include "TChain.h"
#include "TCut.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPRegexp.h"
#include "TROOT.h"

#include <algorithm>
#include <vector>

float GetIntLumi(float lumi, int brun, int bls, int erun, int els)
{
    TChain *c = new TChain("tree");
    //c->Add("/tas05/disk00/jribnik/hunt/dilep_baby/*.root");
    c->Add("/nfs-3/userdata/yanjuntu/hunt/EG_Run2010A-Sep17ReReco_v2_RECO/dilep_baby/*.root");
    c->Add("/nfs-3/userdata/yanjuntu/hunt/Electron_Run2010B-PromptReco-v2_RECO/dilep_baby/*.root");
    //    c->Add("/nfs-3/userdata/yanjuntu/hunt/Mu_Run2010A-Sep17ReReco_v2_RECO/dilep_baby/*.root");
    //c->Add("/nfs-3/userdata/yanjuntu/hunt/Mu_Run2010B-PromptReco-v2_RECO/dilep_baby/*.root");

    TCut c_goodrun(Form("((run>%i&&run<%i)||(run==%i&&ls>=%i)||(run==%i&&ls<=%i))&&goodrun_json(run,ls)", brun, erun, brun, bls, erun, els));
    // goodrun plus events beyond range of goodrun
    // which are not goodrun penalized
    TCut c_goodrunplus(Form("(((run>%i&&run<%i)||(run==%i&&ls>=%i)||(run==%i&&ls<=%i))&&goodrun_json(run,ls))||(run>%i||(run==%i&&ls>%i))", brun, erun, brun, bls, erun, els, erun, erun, els));

    // this is of course particular to dilep babies
    TCut c_notduplicate("! is_duplicate(run,evt,ls,pt1,pt2)");

    // brun:bls -> erun:els
    reset_babydorkidentifier();
    int n_goodrun = c->GetEntries(c_goodrun+c_notduplicate+inclusivez_dilep); 
    // total
    reset_babydorkidentifier();
    int n_total   = c->GetEntries(c_goodrunplus+c_notduplicate+inclusivez_dilep);
    // that which is new
    int n_new     = n_total-n_goodrun;

    float newlumi = ((float)(n_new*lumi))/(float)n_goodrun;
    return lumi+newlumi;
}

float GetIntLumi(float lumi)
{
    int brun = min_run();
    int bls  = min_run_min_lumi();
    int erun = max_run();
    int els  = max_run_max_lumi();

    return GetIntLumi(lumi, brun, bls, erun, els);
}

TH1F* Plot(const char *pfx, TChain *chain, const char *field, TCut sel, TCut presel, float intlumifb, float kfactor,
        unsigned int nbins, float xlo, float xhi, bool integrated, bool isdata)
{
    TCut scale;
    if (! isdata)
        scale = Form("scale1fb*%f*%f", intlumifb, kfactor);
    else
        scale = Form("%f*%f", intlumifb, kfactor);

    char *name = 0;
    if (integrated) name = Form("%s_%s_%s_int", pfx, sel.GetName(), field);
    else name = Form("%s_%s_%s", pfx, sel.GetName(), field);

    char *title = Form("%s_%s_%s, ~%.2f/pb", pfx, sel.GetName(), field, 1e3*intlumifb);

    TH1F *h = 0;
    if (! (h = (TH1F*)gROOT->FindObjectAny(name)))
        h = new TH1F(name, title, nbins, xlo, xhi);

    int brun = min_run();
    int bls  = min_run_min_lumi();
    int erun = max_run();
    int els  = max_run_max_lumi();

    // Used for data only
    TCut c_goodrunplus(Form("(((run>%i&&run<%i)||(run==%i&&ls>=%i)||(run==%i&&ls<=%i))&&goodrun_json(run,ls))||(run>%i||(run==%i&&ls>%i))", brun, erun, brun, bls, erun, els, erun, erun, els));
    TCut c_notduplicate = "";

    // is this an emu, dilep or trilep baby?
    if (chain->GetBranch("pt3")) // trilep
        c_notduplicate = ("! is_duplicate(run,evt,ls,pt1,pt2,pt3)");
    else if (chain->GetBranch("pt2")) // dilep
        c_notduplicate = ("! is_duplicate(run,evt,ls,pt1,pt2)");
    else
        c_notduplicate = ("! is_duplicate(run,evt,ls,pt1)");

    TCut c_presel = "";
    // apply goodrun selection to data where
    // applicable, i.e. the range covered by
    // the goodrun json file
    // also apply BabyDork duplicate removal
    if (isdata)
        c_presel += c_goodrunplus+c_notduplicate+presel;
    else
        c_presel += presel;

    char *draw = Form("%s>>+%s", field, name);
    TCut cut = scale*(c_presel+sel);
    chain->Draw(draw, cut, "goff");

    // Move overflow to the last bin
    float overflow = h->GetBinContent(nbins+1);
    h->SetBinContent(nbins,overflow);
    h->SetBinContent(nbins+1,0.);

    return h;
}

TH1F* Plot(const char *field, TCut sel, float intlumifb, unsigned int nbins, float xlo, float xhi, bool integrated,
        BabySample *bs)
{
    return Plot(bs->pfx(),bs->chain(),field,sel,bs->presel(),intlumifb,bs->kfactor(),nbins,xlo,xhi,integrated,bs->isdata());
}

//
// These should only be used with data, where intlumifb need not be specified and 
// most likely you aren't scaling
//

TH1F* Plot(const char *pfx, TChain *chain, const char *field, TCut sel, TCut presel, float kfactor,
        unsigned int nbins, float xlo, float xhi, bool integrated, bool isdata)
{
    return Plot(pfx,chain,field,sel,presel,1,kfactor,nbins,xlo,xhi,integrated,isdata);
}

TH1F* Plot(const char *pfx, TChain *chain, const char *field, TCut sel, TCut presel,
        unsigned int nbins, float xlo, float xhi, bool integrated, bool isdata)
{
    return Plot(pfx,chain,field,sel,presel,1,1,nbins,xlo,xhi,integrated,isdata);
}

TH1F* Plot(const char *field, TCut sel, unsigned int nbins, float xlo, float xhi, bool integrated,
        BabySample *bs)
{
    return Plot(bs->pfx(),bs->chain(),field,sel,bs->presel(),1.,bs->kfactor(),nbins,xlo,xhi,integrated,bs->isdata());
}

bool sortHistsByIntegral(TH1* h1, TH1* h2)
{
    return h1->Integral() > h2->Integral();
}

// integrate from bin to infinity
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

TCanvas* DrawAll(const char *field, TCut sel, float intlumifb, unsigned int nbins, float xlo, float xhi, bool integrated,
        std::vector<BabySample*> bss)
{
    std::vector<TH1F*> hmcs;
    std::vector<TH1F*> hdatas;
    TH1F* buffer;

    // reset BabyDorkIdentifier here so that
    // it is used for all data BabySamples
    reset_babydorkidentifier();

    for(unsigned int i = 0; i < bss.size(); ++i) {
        if (bss[i]) {
            if (bss[i]->isdata()) {
                buffer = Plot(field,sel,nbins,xlo,xhi,integrated,bss[i]);

                std::vector<TH1F*>::const_iterator it;
                it = find(hdatas.begin(),hdatas.end(),buffer);

                if (it == hdatas.end()) {
                    buffer->SetMarkerColor(bss[i]->color());
                    buffer->SetMarkerStyle(bss[i]->style());
                    hdatas.push_back(buffer);
                }
            } else {
                buffer = Plot(field,sel,intlumifb,nbins,xlo,xhi,integrated,bss[i]);

                std::vector<TH1F*>::const_iterator it;
                it = find(hmcs.begin(),hmcs.end(),buffer);

                if (it == hmcs.end()) {
                    buffer->SetFillColor(bss[i]->color());
                    buffer->SetFillStyle(bss[i]->style());
                    hmcs.push_back(buffer);
                }
            }
        }
    }

    // sort histograms by their total contributions
    // makes for prettier stack plots
    sort(hmcs.begin(), hmcs.end(), sortHistsByIntegral);
    sort(hdatas.begin(), hdatas.end(), sortHistsByIntegral);

    for(unsigned int hmcsIndex = 0; hmcsIndex < hmcs.size(); ++hmcsIndex) {
        for(unsigned int hmcsIndex2 = hmcsIndex+1; hmcsIndex2 < hmcs.size(); ++hmcsIndex2) {
            (hmcs.at(hmcsIndex))->Add(hmcs.at(hmcsIndex2),1);
        }
    }
    for(unsigned int hdatasIndex = 0; hdatasIndex < hdatas.size(); ++hdatasIndex) {
        for(unsigned int hdatasIndex2 = hdatasIndex+1; hdatasIndex2 < hdatas.size(); ++hdatasIndex2) {
            (hdatas.at(hdatasIndex))->Add(hdatas.at(hdatasIndex2),1);
        }
    }

    TH1F *hmc = hmcs[0], *hdata = hdatas[0];
    if (integrated) {
        hmc   = slideIntegrated( hmcs[0]   );
        hdata = slideIntegrated( hdatas[0] );
    }

    char *c1name = 0;
    if (integrated) c1name = Form("c1_%s_%s_int", sel.GetName(), field);
    else c1name = Form("c1_%s_%s", sel.GetName(), field);
    TCanvas *c1 = new TCanvas(c1name);
    c1->Draw();

    hmc->SetTitle(Form("%s, ~%.2f/pb", sel.GetName(), 1e3*intlumifb));
    float ymax = hdata->GetMaximum() > hmc->GetMaximum() ? hdata->GetMaximum() : hmc->GetMaximum();
    hmc->SetMaximum(ymax*1.25);
    hmc->Draw("hist");
    if (integrated)
        hmc->GetXaxis()->SetTitle(Form("integrated %s",field));
    else
        hmc->GetXaxis()->SetTitle(field);

    TLegend* leg = new TLegend(0.82,0.61,0.96,0.93);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetShadowColor(0);

    // For extracting sample prefix
    TPRegexp preg("^([^_]+)_.*$");
    TString  s_pfx("");

    // Add data to legend first so
    // that it's at the top
    hdata->SetMarkerStyle(20);
    if (preg.MatchB(TString(hdata->GetName())))
        s_pfx = ((TObjString*)(preg.MatchS(TString(hdata->GetName()))->At(1)))->GetString();
    leg->AddEntry(hdata, s_pfx.Data(), "lep");

    // And back to MC
    s_pfx = "";
    if (preg.MatchB(TString(hmc->GetName())))
        s_pfx = ((TObjString*)(preg.MatchS(TString(hmc->GetName()))->At(1)))->GetString();
    leg->AddEntry(hmc, s_pfx.Data(), "f");

    for(unsigned int mchisto = 1; mchisto < hmcs.size(); ++mchisto) { // need 1 here, as we only want to plot 1 and more
        if (integrated)
            hmcs.at(mchisto) = slideIntegrated( hmcs.at(mchisto) );
        hmcs.at(mchisto)->Draw("histsame");

        s_pfx = "";
        if (preg.MatchB(TString(hmc->GetName())))
            s_pfx   = ((TObjString*)(preg.MatchS(TString(hmcs.at(mchisto)->GetName()))->At(1)))->GetString();
        leg->AddEntry(hmcs.at(mchisto), s_pfx.Data(), "f");
    }

    // Finally the data
    hdata->Draw("pesame");

    leg->Draw();
    c1->RedrawAxis();

    return c1;
}

TCanvas* DrawAll(const char *field, TCut sel, float intlumifb, unsigned int nbins, float xlo, float xhi, bool integrated,
        BabySample* bs1,  BabySample* bs2,  BabySample* bs3,  BabySample* bs4,  BabySample* bs5,
        BabySample* bs6,  BabySample* bs7,  BabySample* bs8,  BabySample* bs9,  BabySample* bs10,
        BabySample* bs11, BabySample* bs12, BabySample* bs13, BabySample* bs14, BabySample* bs15)
{
    std::vector<BabySample*> tmp;
    BabySample* bss[15] = {bs1,bs2,bs3,bs4,bs5,bs6,bs7,bs8,bs9,bs10,bs11,bs12,bs13,bs14,bs15};
    for(unsigned int i = 0; i < 15; ++i) {
        if (bss[i])
            tmp.push_back(bss[i]);
    }

    return DrawAll(field,sel,intlumifb,nbins,xlo,xhi,integrated,tmp);
}

// Predefines what are most likley the only BabySamples
// one needs for gathering
TCanvas* DrawAll(const char *field, TCut sel, float intlumifb, unsigned int nbins, float xlo, float xhi, bool integrated)
{
    //
    // data babies
    //

    //static BabySample *bs_data_emu    = new BabySample("data","/tas05/disk00/jribnik/hunt/emu_baby/*.root","",1.,true);
    //static BabySample *bs_data_dilep  = new BabySample("data","/tas05/disk00/jribnik/hunt/dilep_baby/*.root","",1.,true);
  static BabySample *bs_data_dilep_1  = new BabySample("data","/nfs-3/userdata/yanjuntu/hunt/EG_Run2010A-Sep17ReReco_v2_RECO/dilep_baby/*.root","",1.,true);
  static BabySample *bs_data_dilep_2  = new BabySample("data","/nfs-3/userdata/yanjuntu/hunt/Electron_Run2010B-PromptReco-v2_RECO/dilep_baby/*.root","",1.,true);
  static BabySample *bs_data_dilep_3  = new BabySample("data","/nfs-3/userdata/yanjuntu/hunt/Mu_Run2010B-PromptReco-v2_RECO/dilep_baby/*.root","",1.,true);
  static BabySample *bs_data_dilep_4  = new BabySample("data","/nfs-3/userdata/yanjuntu/hunt/Mu_Run2010A-Sep17ReReco_v2_RECO/dilep_baby/*.root","",1.,true);
    //static BabySample *bs_data_trilep = new BabySample("data","/tas05/disk00/jribnik/hunt/trilep_baby/*.root","",1.,true);

    //
    // mc babies
    //

    // kfactors
    float kttbarjets = 157.5/165.;
    float ksingletop = 1.;
    float kvvjets    = 1.;
    float kwjets     = 31314./28049.;
    float kzjets     = 3048./2400.;
    float kzll       = 1666./1300.;
    float kdyll      = 3457./2659.;

    // dilep
    static BabySample *bs_ttbarjets_dilep = new BabySample("ttbar","/nfs-3/userdata/yanjuntu/huntmc/TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root","",kttbarjets,false,kRed+1,1001);
    static BabySample *bs_singletop_dilep = new BabySample("tW","/nfs-3/userdata/yanjuntu/huntmc/SingleTop_tWChannel-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root","",ksingletop,false,kMagenta,1001);

    static BabySample *bs_vvjets_dilep    = new BabySample("vvjets","/nfs-3/userdata/yanjuntu/huntmc/VVJets-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root","",kvvjets,false,10,1001);
    static BabySample *bs_wjets_dilep     = new BabySample("wjets","/nfs-3/userdata/yanjuntu/huntmc/WJets-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root","",kwjets,false,kGreen-3,1001);
    static BabySample *bs_ztautau_dilep   = new BabySample("ztautau","/nfs-3/userdata/yanjuntu/huntmc/Ztautau_Spring10-START3X_V26_S09-v1/dilep_baby/*.root","",kzll,false,kAzure+8,1001);

    // Note that a common prefix means
    // a common histogram when used in
    // the same DrawAll, i.e. the five
    // samples below are combined
    static BabySample *bs_zjets_dilep     = new BabySample("zll","/nfs-3/userdata/yanjuntu/huntmc/ZJets-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root","",kzjets,false,kAzure-2,1001);
    static BabySample *bs_zee_dilep       = new BabySample("zll","/nfs-3/userdata/yanjuntu/huntmc/Zee_Spring10-START3X_V26_S09-v1/dilep_baby/*.root","mass<50",kzll,false,kAzure-2,1001);
    static BabySample *bs_zmumu_dilep     = new BabySample("zll","/nfs-3/userdata/yanjuntu/huntmc/Zmumu_Spring10-START3X_V26_S09-v1/dilep_baby/*.root","mass<50",kzll,false,kAzure-2,1001);
    static BabySample *bs_dyee_dilep      = new BabySample("zll","/nfs-3/userdata/yanjuntu/huntmc/DYee_M10to20_Spring10-START3X_V26_S09-v1/dilep_baby/*.root","",kdyll,false,kAzure-2,1001);
    static BabySample *bs_dymumu_dilep    = new BabySample("zll","/nfs-3/userdata/yanjuntu/huntmc/DYmumu_M10to20_Spring10-START3X_V26_S09-v1/dilep_baby/*.root","",kdyll,false,kAzure-2,1001);

    return DrawAll(field,sel,intlumifb,nbins,xlo,xhi,integrated,bs_data_dilep_1,bs_data_dilep_2,bs_data_dilep_3,bs_data_dilep_4,bs_ttbarjets_dilep,bs_singletop_dilep,bs_dyee_dilep,bs_dymumu_dilep,bs_vvjets_dilep,bs_wjets_dilep,bs_zjets_dilep,bs_zee_dilep,bs_zmumu_dilep,bs_ztautau_dilep);
}

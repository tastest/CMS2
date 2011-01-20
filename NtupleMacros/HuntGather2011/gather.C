#include "BabyDorkIdentifier.h"
#include "BabySample.h"
#include "cuts.h"
#include "gather.h"
#include "../Tools/goodrun.h"

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

float GetIntLumi(float lumi, int brun, int bls, int erun, int els)
{
    TChain *c = new TChain("tree");
    //c->Add("/tas05/disk00/jribnik/hunt/dilep_baby/*.root");
    c->Add("/nfs-3/userdata/yanjuntu/hunt/EG_Run2010A-Sep17ReReco_v2_RECO/V03-06-14/dilep_baby/*.root");
    c->Add("/nfs-3/userdata/yanjuntu/hunt/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14-00/dilep_baby/*.root");
    c->Add("/nfs-3/userdata/yanjuntu/hunt/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14/dilep_baby/*.root");
    c->Add("/nfs-3/userdata/yanjuntu/hunt/Mu_Run2010A-Sep17ReReco_v2_RECO/V03-06-14/dilep_baby/*.root");
    c->Add("/nfs-3/userdata/yanjuntu/hunt/Mu_Run2010B-PromptReco-v2_RECO/V03-06-14-00/dilep_baby/*.root");
    c->Add("/nfs-3/userdata/yanjuntu/hunt/Mu_Run2010B-PromptReco-v2_RECO/V03-06-14/dilep_baby/*.root");

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

    std::cout << "in the good run list for " << lumi << " there are " << n_goodrun << std::endl;
    std::cout << "out of the good run list there are " << n_new << std::endl;
    std::cout << "that means a new lumi of " << newlumi << std::endl;


    return lumi+newlumi;
}

float GetIntLumi(float lumi)
{
    int brun = min_run();
    int bls  = min_run_min_lumi();
    int erun = max_run();
    int els  = max_run_max_lumi();
    std::cout << brun << ", " << bls << " : " << erun << ", " << els << std::endl;

    return GetIntLumi(lumi, brun, bls, erun, els);
}

TH1F* Plot(const char *pfx, const char *pfx2,TChain *chain, TCut field, TCut sel, TCut presel, float intlumifb, float kfactor,
        unsigned int nbins, float xlo, float xhi, bool integrated, bool isdata, unsigned int isfx)
{
    TCut scale;
    if (! isdata)
        scale = Form("scale1fb*%f*%f", intlumifb, kfactor);
    else
        scale = Form("%f*%f", intlumifb, kfactor);

    char *name = 0;
    if (integrated) name = Form("%s_%s_%s_%i_int", pfx, sel.GetName(), field.GetName(), isfx);
    else name = Form("%s_%s_%s_%i", pfx, sel.GetName(), field.GetName(), isfx);

    char *title = Form("%s_%s_%s, ~%.2f/pb", pfx, sel.GetName(), field.GetName(), 1e3*intlumifb);

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

    char *draw = Form("%s>>+%s", field.GetTitle(), name);
    // TCut cut = scale*(c_presel+sel);
    TCut cut = c_presel+sel;

    // Okay, time for someting sneaky. At the
    // time of this commit I anticipate YJ
    // adding a chain->Scan immediately after
    // the chain->Draw. She will do this only
    // for data. The Draw and Scan use exactly
    // the same selection, so let's first
    // create a TEventList with the selection
    // and use that for both. Be sure to go
    // back to the original TEventList when
    // you're done.
    TEventList* orig_elist = (TEventList*)chain->GetEventList()->Clone();
    chain->Draw(">>evtlist", cut, "goff");
    TEventList* elist = (TEventList*)gDirectory->Get("evtlist");
    chain->SetEventList(elist);
    chain->Draw(draw, scale, "goff");
    // chain->Scan(var1:var2...", "no cuts")
    // All done now, revert to original event list
   
    if (isdata){
      //reset_babydorkidentifier();                                                                                                                                              
      TTreePlayer *tp = (TTreePlayer*)chain->GetPlayer();
      tp->SetScanRedirect(kTRUE);
      //tp->SetScanFileName(Form("%s_%s_%s_%i.out", pfx, sel.GetName(), field.GetName(), isfx));                                                                                 
      tp->SetScanFileName(Form("%s_%s_%i_%s.out", sel.GetName(), field.GetName(), isfx,pfx2));
      if (chain->GetBranch("pt3")) // trilep                                                                                                                                     
        chain->Scan("dataset:run:ls:evt:pt1:pt2:pt3", "", "colsize=50");
      else if (chain->GetBranch("pt2")) // dilep                                                                                                                                 
        chain->Scan("dataset:run:ls:evt:pt1:pt2", "", "colsize=50");
      else
        chain->Scan("dataset:run:ls:evt:pt1", "", "colsize=50");
    }
    chain->SetEventList(orig_elist);


    // Move overflow to the last bin
    float overflow = h->GetBinContent(nbins+1);
    h->SetBinContent(nbins,overflow);
    h->SetBinContent(nbins+1,0.);

    return h;
}

TH1F* Plot(TCut field, TCut sel, TCut presel, float intlumifb, unsigned int nbins, float xlo, float xhi, bool integrated,
        BabySample *bs, unsigned int isfx)
{
  return Plot(bs->pfx(),bs->pfx2(),bs->chain(),field,sel,presel+bs->presel(),intlumifb,bs->kfactor(),nbins,xlo,xhi,integrated,bs->isdata(),isfx);
}

//
// These should only be used with data, where intlumifb need not be specified and 
// most likely you aren't scaling
//

TH1F* Plot(const char *pfx, const char *pfx2,TChain *chain, TCut field, TCut sel, TCut presel, float kfactor,
        unsigned int nbins, float xlo, float xhi, bool integrated, bool isdata, unsigned int isfx)
{
  return Plot(pfx,pfx2,chain,field,sel,presel,1,kfactor,nbins,xlo,xhi,integrated,isdata,isfx);
}

TH1F* Plot(const char *pfx, const char *pfx2,TChain *chain, TCut field, TCut sel, TCut presel,
        unsigned int nbins, float xlo, float xhi, bool integrated, bool isdata, unsigned int isfx)
{
  return Plot(pfx,pfx2,chain,field,sel,presel,1,1,nbins,xlo,xhi,integrated,isdata,isfx);
}

TH1F* Plot(TCut field, TCut sel, TCut presel, unsigned int nbins, float xlo, float xhi, bool integrated,
        BabySample *bs, unsigned int isfx)
{
  // return Plot(bs->pfx(),bs->chain(),field,sel,presel+bs->presel(),1.,bs->kfactor(),nbins,xlo,xhi,integrated,bs->isdata(),isfx);
  return Plot(bs->pfx(),bs->pfx2(),bs->chain(),field,sel,presel+bs->presel(),1.,bs->kfactor(),nbins,xlo,xhi,integrated,bs->isdata(),isfx);
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

TCanvas* DrawAll(TCut field, const char *savename, TCut sel, TCut presel, float intlumifb, unsigned int nbins, float xlo, float xhi, bool integrated,
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
                buffer = Plot(field,sel,presel,nbins,xlo,xhi,integrated,bss[i],gDrawAllCount);

                std::vector<TH1F*>::const_iterator it;
                it = find(hdatas.begin(),hdatas.end(),buffer);

                if (it == hdatas.end()) {
                    buffer->SetMarkerColor(bss[i]->color());
                    buffer->SetMarkerStyle(bss[i]->style());
                    hdatas.push_back(buffer);
                }
            } else {
                buffer = Plot(field,sel,presel,intlumifb,nbins,xlo,xhi,integrated,bss[i],gDrawAllCount);

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

    TCanvas *c1 = new TCanvas(savename);
    c1->Draw();

    if (!strcmp("CUT", field.GetName())) 
        field.SetName(field.GetTitle());

    hmc->SetTitle(Form("%s, ~%.2f/pb", sel.GetName(), 1e3*intlumifb));
    float ymax = hdata->GetMaximum() > hmc->GetMaximum() ? hdata->GetMaximum()+2*sqrt(hdata->GetMaximum()) : hmc->GetMaximum()*1.25;
    hmc->SetMaximum(ymax);
    hmc->Draw("hist");
    if (integrated)
        hmc->GetXaxis()->SetTitle(Form("integrated %s",field.GetName()));
    else
        hmc->GetXaxis()->SetTitle(field.GetName());

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

    gDrawAllCount++;
    return c1;
}

TCanvas* DrawAll(TCut field, const char *savename, TCut sel, TCut presel, float intlumifb, unsigned int nbins, float xlo, float xhi, bool integrated,
		 BabySample* bs1,  BabySample* bs2,  BabySample* bs3,  BabySample* bs4,  BabySample* bs5,
		 BabySample* bs6,  BabySample* bs7,  BabySample* bs8,  BabySample* bs9,  BabySample* bs10,
		 BabySample* bs11, BabySample* bs12, BabySample* bs13, BabySample* bs14, BabySample* bs15,
		 BabySample* bs16, BabySample* bs17, BabySample* bs18, BabySample* bs19, BabySample* bs20)
{
    std::vector<BabySample*> tmp;
    BabySample* bss[20] = {bs1,bs2,bs3,bs4,bs5,bs6,bs7,bs8,bs9,bs10,bs11,bs12,bs13,bs14,bs15,bs16,bs17,bs18,bs19,bs20};
    for(unsigned int i = 0; i < 20; ++i) {
        if (bss[i])
            tmp.push_back(bss[i]);
    }

    return DrawAll(field,savename,sel,presel,intlumifb,nbins,xlo,xhi,integrated,tmp);
}

// Predefines what are most likley the only BabySamples
// one needs for gathering
TCanvas* DrawAll(TCut field, const char *savename, TCut sel, TCut presel, float intlumifb, unsigned int nbins, float xlo, float xhi, bool integrated)
{
    // apply base_dilep selection to dilep babies
    // this speeds things up

    //
    // data babies
    //

  static BabySample *bs_data_dilep_1  = new BabySample("data","Electron1","/nfs-3/userdata/yanjuntu/hunt/EG_Run2010A-Sep17ReReco_v2_RECO/V03-06-14/dilep_baby/*.root",base_dilep,1.,true);
  static BabySample *bs_data_dilep_2  = new BabySample("data","Electron2","/nfs-3/userdata/yanjuntu/hunt/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14-00/dilep_baby/*.root",base_dilep,1.,true);
  static BabySample *bs_data_dilep_3  = new BabySample("data","Electron3","/nfs-3/userdata/yanjuntu/hunt/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14/dilep_baby/*.root",base_dilep,1.,true);
  static BabySample *bs_data_dilep_4  = new BabySample("data","Mu1","/nfs-3/userdata/yanjuntu/hunt/Mu_Run2010A-Sep17ReReco_v2_RECO/V03-06-14/dilep_baby/*.root",base_dilep,1.,true); 
  static BabySample *bs_data_dilep_5  = new BabySample("data","Mu2","/nfs-3/userdata/yanjuntu/hunt/Mu_Run2010B-PromptReco-v2_RECO/V03-06-14-00/dilep_baby/*.root",base_dilep,1.,true);
  static BabySample *bs_data_dilep_6  = new BabySample("data","Mu3","/nfs-3/userdata/yanjuntu/hunt/Mu_Run2010B-PromptReco-v2_RECO/V03-06-14/dilep_baby/*.root",base_dilep,1.,true);

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

    // stitching
    TCut cut_notau("ngentaus==0");
    TCut cut_tau("ngentaus==2");
    TCut cut_stitch_mass("mass<50");
    TCut cut_dilep_stitch_mass = base_dilep + cut_stitch_mass;
    TCut cut_dilep_notau = base_dilep + cut_notau;
    TCut cut_dilep_tau = base_dilep + cut_tau;

    // dilep
    static BabySample *bs_ttbarjets_dilep = new BabySample("ttbar","mc","/nfs-3/userdata/yanjuntu/huntmc/TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root",base_dilep,kttbarjets,false,kRed+1,1001);
    static BabySample *bs_singletop_dilep = new BabySample("tW","mc","/nfs-3/userdata/yanjuntu/huntmc/SingleTop_tWChannel-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root",base_dilep,ksingletop,false,kMagenta,1001);

    static BabySample *bs_vvjets_dilep    = new BabySample("vvjets","mc","/nfs-3/userdata/yanjuntu/huntmc/VVJets-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root",base_dilep,kvvjets,false,10,1001);
    static BabySample *bs_wjets_dilep     = new BabySample("wjets","mc","/nfs-3/userdata/yanjuntu/huntmc/WJets-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root",base_dilep,kwjets,false,kGreen-3,1001);

    static BabySample *bs_zjetstautau_dilep = new BabySample("ztautau","mc","/nfs-3/userdata/yanjuntu/huntmc/ZJets-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root", cut_dilep_tau,kzll,false,kAzure+8,1001);
    static BabySample *bs_ztautau_dilep     = new BabySample("ztautau","mc","/nfs-3/userdata/yanjuntu/huntmc/Ztautau_Spring10-START3X_V26_S09-v1/dilep_baby/*.root",cut_dilep_stitch_mass,kzll,false,kAzure+8,1001);

    // Note that a common prefix means
    // a common histogram when used in
    // the same DrawAll, i.e. the five
    // samples below are combined

    static BabySample *bs_zjets_dilep     = new BabySample("zll","mc","/nfs-3/userdata/yanjuntu/huntmc/ZJets-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root",cut_dilep_notau,kzjets,false,kAzure-2,1001);
    static BabySample *bs_zee_dilep       = new BabySample("zll","mc","/nfs-3/userdata/yanjuntu/huntmc/Zee_Spring10-START3X_V26_S09-v1/dilep_baby/*.root",cut_dilep_stitch_mass,kzll,false,kAzure-2,1001);
    static BabySample *bs_zmumu_dilep     = new BabySample("zll","mc","/nfs-3/userdata/yanjuntu/huntmc/Zmumu_Spring10-START3X_V26_S09-v1/dilep_baby/*.root",cut_dilep_stitch_mass,kzll,false,kAzure-2,1001);
    static BabySample *bs_dyee_dilep      = new BabySample("zll","mc","/nfs-3/userdata/yanjuntu/huntmc/DYee_M10to20_Spring10-START3X_V26_S09-v1/dilep_baby/*.root",base_dilep,kdyll,false,kAzure-2,1001);
    static BabySample *bs_dymumu_dilep    = new BabySample("zll","mc","/nfs-3/userdata/yanjuntu/huntmc/DYmumu_M10to20_Spring10-START3X_V26_S09-v1/dilep_baby/*.root",base_dilep,kdyll,false,kAzure-2,1001);

    std::vector<BabySample*> babyVector;
    babyVector.push_back(bs_data_dilep_1);
    babyVector.push_back(bs_data_dilep_2);
    babyVector.push_back(bs_data_dilep_3);
    babyVector.push_back(bs_data_dilep_4);
    babyVector.push_back(bs_data_dilep_5);
    babyVector.push_back(bs_data_dilep_6);
    babyVector.push_back(bs_ttbarjets_dilep);
    babyVector.push_back(bs_singletop_dilep);
    babyVector.push_back(bs_dyee_dilep);
    babyVector.push_back(bs_dymumu_dilep);
    babyVector.push_back(bs_vvjets_dilep);
    babyVector.push_back(bs_wjets_dilep);
    babyVector.push_back(bs_zjetstautau_dilep);
    babyVector.push_back(bs_zjets_dilep);
    babyVector.push_back(bs_zee_dilep);
    babyVector.push_back(bs_zmumu_dilep);
    babyVector.push_back(bs_ztautau_dilep);

    return DrawAll(field,savename,sel,presel,intlumifb,nbins,xlo,xhi,integrated,babyVector);
}

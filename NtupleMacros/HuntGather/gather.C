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

float GetIntLumi(float lumi, int brun, int bls, int erun, int els)
{
    TChain *c = new TChain("tree");
    c->Add("/tas05/disk00/jribnik/hunt/dilep_baby/*.root");

    TCut c_goodrun(Form("((run>%i&&run<%i)||(run==%i&&ls>=%i)||(run==%i&&ls<=%i))&&goodrun_json(run,ls)", brun, erun, brun, bls, erun, els));
    // goodrun plus events beyond range of goodrun
    // which are not goodrun penalized
    TCut c_goodrunplus(Form("(((run>%i&&run<%i)||(run==%i&&ls>=%i)||(run==%i&&ls<=%i))&&goodrun_json(run,ls))||(run>%i||(run==%i&&ls>%i))", brun, erun, brun, bls, erun, els, erun, erun, els));

    // brun:bls -> erun:els
    int n_goodrun = c->GetEntries(c_goodrun+inclusivez_dilep); 
    // total
    int n_total   = c->GetEntries(c_goodrunplus+inclusivez_dilep);
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
        unsigned int nbins, float xlo, float xhi, bool usegoodrun)
{
    TCut scale;
    // mc has scale1fb, data does not
    // JAKE is this foolproof? rather than
    // bool usegoodrun should we bool
    // isdata and force application of the
    // goodrun selection?
    if (chain->GetBranch("scale1fb"))
        scale = Form("scale1fb*%f*%f", intlumifb, kfactor);
    else
        scale = Form("%f*%f", intlumifb, kfactor);

    char *name = Form("%s_%s_%s", pfx, sel.GetName(), field);
    TH1F *h = 0;
    if (! (h = (TH1F*)gROOT->FindObjectAny(name)))
        h = new TH1F(name, name, nbins, xlo, xhi);

    int brun = min_run();
    int bls  = min_run_min_lumi();
    int erun = max_run();
    int els  = max_run_max_lumi();

    TCut c_goodrunplus(Form("(((run>%i&&run<%i)||(run==%i&&ls>=%i)||(run==%i&&ls<=%i))&&goodrun_json(run,ls))||(run>%i||(run==%i&&ls>%i))", brun, erun, brun, bls, erun, els, erun, erun, els));

    TCut c_presel = "";
    // apply goodrun selection to data where
    // applicable, i.e. the range covered by
    // the goodrun json file
    if (usegoodrun)
        c_presel += c_goodrunplus+presel;
    else
        c_presel += presel;

    char *draw = Form("%s>>+%s", field, name);
    TCut cut = scale*(c_presel+sel);
    chain->Draw(draw, cut, "goff");

    return h;
}

TH1F* Plot(const char *field, TCut sel, float intlumifb, unsigned int nbins, float xlo, float xhi, bool usegoodrun,
        BabySample *bs)
{
    return Plot(bs->pfx(),bs->chain(),field,sel,bs->presel(),intlumifb,bs->kfactor(),nbins,xlo,xhi,usegoodrun);
}

//
// These should only be used with data, where intlumifb need not be specified and 
// most likely you aren't scaling
//

TH1F* Plot(const char *pfx, TChain *chain, const char *field, TCut sel, TCut presel, float kfactor,
        unsigned int nbins, float xlo, float xhi, bool usegoodrun)
{
    return Plot(pfx,chain,field,sel,presel,1,kfactor,nbins,xlo,xhi,usegoodrun);
}

TH1F* Plot(const char *pfx, TChain *chain, const char *field, TCut sel, TCut presel,
        unsigned int nbins, float xlo, float xhi, bool usegoodrun)
{
    return Plot(pfx,chain,field,sel,presel,1,1,nbins,xlo,xhi,usegoodrun);
}

TH1F* Plot(const char *field, TCut sel, unsigned int nbins, float xlo, float xhi, bool usegoodrun,
        BabySample *bs)
{
    return Plot(bs->pfx(),bs->chain(),field,sel,bs->presel(),1.,bs->kfactor(),nbins,xlo,xhi,usegoodrun);
}

TH1F* slideIntegrated(TH1F* integrateThis)
{
    TString name = integrateThis->GetName();
    name.Append("_int");
    int NbinsX = integrateThis->GetNbinsX();
    TH1F* integrated = (TH1F*)integrateThis->Clone(name);
    integrated->Reset();

    // integrate rate above a certain pt value (bin)
    for(int i = 1; i <= NbinsX; ++i) {
        // don't forget the overflow!
        float integral = integrateThis->Integral(i,NbinsX+1);
        if(integral!=0.) integrated->SetBinContent(i,integral);
        if(integral!=0.) integrated->SetBinError(i,0);
        //        if(integral!=0.) integrated->SetBinError(i,TMath::Sqrt(integral));
    }

    return integrated;
}

TCanvas* DrawAll(const char *field, TCut sel, float intlumifb, unsigned int nbins, float xlo, float xhi, bool integrated,
        BabySample* bs1,  BabySample* bs2,  BabySample* bs3,  BabySample* bs4,  BabySample* bs5,
        BabySample* bs6,  BabySample* bs7,  BabySample* bs8,  BabySample* bs9,  BabySample* bs10,
        BabySample* bs11, BabySample* bs12, BabySample* bs13, BabySample* bs14, BabySample* bs15)
{
    std::vector<TH1F*> hmcs;
    std::vector<TH1F*> hdatas;
    TH1F* buffer;

    BabySample* bss[15] = {bs1,bs2,bs3,bs4,bs5,bs6,bs7,bs8,bs9,bs10,bs11,bs12,bs13,bs14,bs15};
    for(unsigned int i = 0; i < 15; ++i) {
        if (bss[i]) {
            if (bss[i]->isdata()) {
                buffer = Plot(field,sel,nbins,xlo,xhi,1,bss[i]);

                std::vector<TH1F*>::const_iterator it;
                it = find(hdatas.begin(),hdatas.end(),buffer);

                if (it == hdatas.end()) {
                    buffer->SetMarkerColor(bss[i]->color());
                    buffer->SetMarkerStyle(bss[i]->style());
                    hdatas.push_back(buffer);
                }
            } else {
                buffer = Plot(field,sel,intlumifb,nbins,xlo,xhi,0,bss[i]);

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

    TLegend* leg = new TLegend(0.85,0.622881,0.997126,0.98,NULL,"brNDC");
    leg->SetLineColor (1);
    leg->SetLineStyle (1);
    leg->SetLineWidth (1);
    leg->SetFillColor (10);
    leg->SetBorderSize(1);

    // mini Hstack function: (should be replaced with something more adequate)
    for(unsigned int mchisto = 0; mchisto < hmcs.size(); ++mchisto) {
        for(unsigned int mchisto2 = mchisto+1; mchisto2 < hmcs.size(); ++mchisto2) {
            (hmcs.at(mchisto))->Add(hmcs.at(mchisto2),1);
        }
    }
    for(unsigned int datahisto = 0; datahisto < hdatas.size(); ++datahisto) {
        for(unsigned int datahisto2 = datahisto+1; datahisto2 < hdatas.size(); ++datahisto2) {
            (hdatas.at(datahisto))->Add(hdatas.at(datahisto2),1);
        }
    }

    TH1F *hmc = hmcs[0], *hdata = hdatas[0];
    if (integrated) {
        hmc   = slideIntegrated( hmcs[0]   );
        hdata = slideIntegrated( hdatas[0] );
    }

    TPRegexp preg("^([^_]+)_([^_]+_[^_]+)_([^_]+)(_int)?$");
    TString  s_pfx(""),s_sel(""),s_field("");
    if (preg.MatchB(TString(hmc->GetName()))) {
        s_sel   = ((TObjString*)(preg.MatchS(TString(hmc->GetName()))->At(2)))->GetString();
        s_field = ((TObjString*)(preg.MatchS(TString(hmc->GetName()))->At(3)))->GetString();
    }

    hmc->SetTitle(s_sel.Data());
    float ymax = hdata->GetMaximum() > hmc->GetMaximum() ? hdata->GetMaximum() : hmc->GetMaximum();
    hmc->SetMaximum(ymax*1.1);

    TCanvas *c1 = new TCanvas("c1");
    c1->Draw();
    // Draw the initial MC histo
    hmc->Draw();
    hmc->GetXaxis()->SetTitle(s_field.Data());

    // But add data to legend first
    // so that it's at the top
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

        hmcs.at(mchisto)->Draw("same");

        s_pfx = "";
        if (preg.MatchB(TString(hmc->GetName())))
            s_pfx   = ((TObjString*)(preg.MatchS(TString(hmcs.at(mchisto)->GetName()))->At(1)))->GetString();

        leg->AddEntry(hmcs.at(mchisto), s_pfx.Data(), "f");
    }

    hdata->Draw("pesames");

    leg->Draw();
    c1->RedrawAxis();

    return c1;
}

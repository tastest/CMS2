#include "cuts.h"
#include "gather.h"
#include "goodrun.h"

#include "TChain.h"
#include "TCut.h"
#include "TH1F.h"

float GetIntLumi(float lumi, int brun, int bls, int erun, int els)
{
    TChain *c = new TChain("tree");
    c->Add("/tas05/disk00/jribnik/hunt/dilep_baby/*.root");

    TCut c_goodrun("goodrun_json(run,ls)");
    TCut c_runls(Form("(run>%i&&run<%i)||(run==%i&&ls>=%i)||(run==%i&&ls<=%i)", brun, erun, brun, bls, erun, els));

    // brun:bls -> erun:els
    int n_runls = c->GetEntries(c_goodrun+c_runls+inclusivez_dilep); 
    // total
    int n_total = c->GetEntries(c_goodrun+inclusivez_dilep);
    // that which is new
    int n_new   = n_total-n_runls;

    float newlumi = ((float)(n_new*lumi))/(float)n_runls;
    return lumi+newlumi;
}

float GetIntLumi(float lumi)
{
    TChain *c = new TChain("tree");
    c->Add("/tas05/disk00/jribnik/hunt/dilep_baby/*.root");

    int brun = min_run();
    int bls  = min_run_min_lumi();
    int erun = max_run();
    int els  = max_run_max_lumi();

    TCut c_goodrun("goodrun_json(run,ls)");
    TCut c_runls(Form("(run>%i&&run<%i)||(run==%i&&ls>=%i)||(run==%i&&ls<=%i)", brun, erun, brun, bls, erun, els));

    // brun:bls -> erun:els
    int n_runls = c->GetEntries(c_goodrun+c_runls+inclusivez_dilep); 
    // total
    int n_total = c->GetEntries(c_goodrun+inclusivez_dilep);
    // that which is new
    int n_new   = n_total-n_runls;

    float newlumi = ((float)(n_new*lumi))/(float)n_runls;
    return lumi+newlumi;
}

TH1F* Plot(const char *pfx, TChain *chain, const char *field, TCut sel, TCut presel, float intlumiinfb, float kfactor,
        unsigned int nbins, float xlo, float xhi)
{
    TCut scale;
    // mc has scale1fb, data does not
    if (chain->GetBranch("scale1fb"))
        scale = Form("scale1fb*%f*%f", intlumiinfb, kfactor);
    else
        scale = Form("%f*%f", intlumiinfb, kfactor);

    char *name = Form("%s_%s_%s", pfx, sel.GetName(), field);
    TH1F *h = new TH1F(name, name, nbins, xlo, xhi);

    char *draw = Form("%s>>+%s", field, name);
    TCut cut = scale*(presel+sel);
    chain->Draw(draw, cut, "goff");

    return h;
}

//
// These should only be used with data, where intlumiinfb need not be specified and 
// most likely you aren't scaling
//

TH1F* Plot(const char *pfx, TChain *chain, const char *field, TCut sel, TCut presel, float kfactor,
        unsigned int nbins, float xlo, float xhi)
{
    return Plot(pfx,chain,field,sel,presel,1,kfactor,nbins,xlo,xhi);
}

TH1F* Plot(const char *pfx, TChain *chain, const char *field, TCut sel, TCut presel,
        unsigned int nbins, float xlo, float xhi)
{
    return Plot(pfx,chain,field,sel,presel,1,1,nbins,xlo,xhi);
}

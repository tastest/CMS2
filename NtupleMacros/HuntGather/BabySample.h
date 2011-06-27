#ifndef BABYSAMPLE_H
#define BABYSAMPLE_H

#include "TChain.h"
#include "TCut.h"
#include "TDirectory.h"
#include "TEventList.h"

#include <cstdio>

class BabySample
{
    public:
        BabySample() {}
        //BabySample(const char *pfx, const char* babies, TCut presel, float kfactor, bool isdata, Color_t color = kBlack, Style_t style = 20)
	BabySample(const char *pfx,const char *pfx2, const char* babies, TCut presel, float kfactor, bool isdata, Color_t color = kBlack, Style_t style = 20)
        {
            color_ = color;
            style_ = style;
            chain_ = new TChain("tree");
            chain_->Add(babies);
            presel_ = presel;
            pfx_ = pfx;
	    pfx2_ = pfx2;
            isdata_ = isdata;
            kfactor_ = kfactor;

            if (strcmp(presel.GetTitle(),"")) {
                chain_->Draw(">>elist", presel);
                elist_ = (TEventList*)gDirectory->Get("elist")->Clone();
                chain_->SetEventList(elist_);
            }
        }
        ~BabySample() {}

        Color_t color()   const { return color_; }
        Style_t style()   const { return style_; }
        TChain* chain()   const { return chain_; }
        TCut    presel()  const { return presel_;}
        const char* pfx() const { return pfx_.Data(); }
	const char* pfx2() const { return pfx2_.Data(); }
        bool    isdata()  const { return isdata_;  }
        float   kfactor() const { return kfactor_; }

    private:
        Color_t color_;
        Style_t style_;
        TChain* chain_;
        TCut    presel_;
        TEventList* elist_;
        TString pfx_;
	TString pfx2_;

        bool  isdata_;
        float kfactor_;
};

#endif

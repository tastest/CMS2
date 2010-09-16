#ifndef BABYSAMPLE_H
#define BABYSAMPLE_H

#include "TChain.h"
#include "TCut.h"

class BabySample
{
    public:
        BabySample() : chain_(0), presel_(""), pfx_(""), isdata_(0), kfactor_(1.), color_(kBlack) {}
        BabySample(char *pfx, TChain* chain, TCut presel, float kfactor, bool isdata, Color_t color = kBlack, Style_t style = 20)
        {
            chain_ = chain;
            presel_ = presel;
            pfx_ = pfx;
            isdata_ = isdata;
            kfactor_ = kfactor;
            color_ = color;
            style_ = style;
        }
        BabySample(char *pfx, char* babies, TCut presel, float kfactor, bool isdata, Color_t color = kBlack, Style_t style = 20)
        {
            chain_ = new TChain("tree");
            chain_->Add(babies);
            presel_ = presel;
            pfx_ = pfx;
            isdata_ = isdata;
            kfactor_ = kfactor;
            color_ = color;
            style_ = style;
        }
        ~BabySample() {}

        Color_t color()   { return color_;      }
        Style_t style()   { return style_;      }
        TChain* chain()   { return chain_;      }
        TCut    presel()  { return presel_;     }
        const char* pfx() { return pfx_.Data(); }
        bool    isdata()  { return isdata_;     }
        float   kfactor() { return kfactor_;    }

    private:
        TChain* chain_;
        TCut    presel_;
        TString pfx_;
        bool  isdata_;
        float kfactor_;
        Color_t color_;
        Style_t style_;
};

#endif

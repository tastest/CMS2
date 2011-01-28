
#ifndef BABYSAMPLE_H
#define BABYSAMPLE_H

#include "TChain.h"
#include "TCut.h"
#include "TDirectory.h"
#include "TEventList.h"

#include <cstdio>
#include "SampleType.h"

class BabySample
{
    public:
        BabySample() {}
        BabySample(const char *pfx,const char *pfx2, const char* babies, TCut presel, float kfactor, SampleType type, Color_t color = kBlack, Style_t style = 20);
        ~BabySample() {}

        void    add(const char* babies);
        Color_t color()   const { return color_; }
        Style_t style()   const { return style_; }
        TChain* chain()   const { return chain_; }
        TCut    presel()  const { return presel_;}
        const char* pfx() const { return pfx_.Data(); }
        const char* pfx2() const { return pfx2_.Data(); }
        SampleType    type()  const { return type_;  }
        float   kfactor() const { return kfactor_; }

    private:
        Color_t color_;
        Style_t style_;
        TChain* chain_;
        TCut    presel_;
        TEventList* elist_;
        TString pfx_;
        TString pfx2_;
        SampleType  type_;
        float kfactor_;
};

#endif


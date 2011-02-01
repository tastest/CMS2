
#include "BabySample.h"
#include <iostream>

BabySample::BabySample(const char *pfx,const char *pfx2, const char* babies, TCut presel, float kfactor, SampleType type, Color_t color, Style_t style)
{
    color_ = color;
    style_ = style;
    chain_ = new TChain("tree");
    chain_->Add(babies);
    presel_ = presel;
    pfx_ = pfx;
    pfx2_ = pfx2;
    type_ = type;
    kfactor_ = kfactor;

// this should be called from elsewhere
//    if (strcmp(presel.GetTitle(),"")) {
//        chain_->Draw(">>elist", presel);
//        elist_ = (TEventList*)gDirectory->Get("elist")->Clone();
//        chain_->SetEventList(elist_);
//    }
}

void BabySample::add(const char* babies) 
{
    std::cout << "adding: " << babies << std::endl;
    chain_->Add(babies);
}


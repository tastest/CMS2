
#include "BabySample.h"

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
    elist_ = 0;

}

void BabySample::add(const char* babies) 
{
    chain_->Add(babies);
}

void BabySample::setEventList(TCut cut)
{

    chain_->SetEventList(0);
    chain_->Draw(">>elist", cut + presel_, "goff");
    elist_ = (TEventList*)gDirectory->Get("elist")->Clone();
    chain_->SetEventList(elist_);

}

void BabySample::resetEventList() 
{
    chain_->SetEventList(0);    
    chain_->SetEventList(elist_);
}


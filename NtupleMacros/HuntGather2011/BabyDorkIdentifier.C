#include "BabyDorkIdentifier.h"

#ifndef __CINT__

#include "TMath.h"
#include <set>

struct BabyDorkIdentifier {
    BabyDorkIdentifier (unsigned int run, unsigned int lumi, unsigned int event, float pt1, float pt2 = -999999., float pt3 = -999999.);
    unsigned int run_, lumi_, event_;
    float pt1_, pt2_, pt3_;
    bool operator < (const BabyDorkIdentifier &) const;
    bool operator == (const BabyDorkIdentifier &) const;
};

BabyDorkIdentifier::BabyDorkIdentifier (unsigned int run, unsigned int lumi, unsigned int event, float pt1, float pt2, float pt3)
     : run_(run), lumi_(lumi), event_(event), pt1_(pt1), pt2_(pt2), pt3_(pt3)
{}

bool BabyDorkIdentifier::operator < (const BabyDorkIdentifier &other) const
{
    if (run_   != other.run_)
        return run_ < other.run_;
    if (lumi_  != other.lumi_)
        return lumi_ < other.lumi_;
    if (event_ != other.event_)
        return event_ < other.event_;
    if (TMath::Abs(pt1_-other.pt1_) > 1e-6*fabs(pt1_))
      return pt1_ < other.pt1_;
    if (TMath::Abs(pt2_-other.pt2_) > 1e-6*fabs(pt2_))
      return pt2_ < other.pt2_;
    if (TMath::Abs(pt3_-other.pt3_) > 1e-6*fabs(pt3_))
      return pt3_ < other.pt3_;
    return false;
}

bool BabyDorkIdentifier::operator == (const BabyDorkIdentifier &other) const
{
    if (run_   != other.run_)
        return false;
    if (lumi_  != other.lumi_)
        return false;
    if (event_ != other.event_)
        return false;
    if (TMath::Abs(pt1_-other.pt1_) > 1e-6*fabs(pt1_))
        return false;
    if (TMath::Abs(pt2_-other.pt2_) > 1e-6*fabs(pt2_))
        return false;
    if (TMath::Abs(pt3_-other.pt3_) > 1e-6*fabs(pt3_))
        return false;
    return true;
}

static std::set<BabyDorkIdentifier> already_seen_;
bool is_duplicate (const BabyDorkIdentifier &id)
{
    std::pair<std::set<BabyDorkIdentifier>::const_iterator, bool> ret =
        already_seen_.insert(id);
    return !ret.second;
}

bool is_duplicate(unsigned int run, unsigned int lumi, unsigned int event, float pt1, float pt2, float pt3)
{
    return is_duplicate(BabyDorkIdentifier(run,lumi,event,pt1,pt2,pt3));
}

void reset_babydorkidentifier()
{
    already_seen_.clear();
}

#endif

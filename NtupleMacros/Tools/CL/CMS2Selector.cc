// $Id: CMS2Selector.cc,v 1.1 2010/04/16 19:01:22 jmuelmen Exp $

#include <stdio.h>
#include "CORE/CMS2.h"
#include "CMS2Selector.h"

CMS2Selector::CMS2Selector ()
     : TSelectorDraw()
{

}

void CMS2Selector::Begin (TTree *t)
{
     cms2.Init(t);
     TSelectorDraw::Begin(t);
}

void CMS2Selector::Init (TTree *t)
{
     cms2.Init(t);
     TSelectorDraw::Init(t);
}

Bool_t CMS2Selector::Notify ()
{
     printf("Notify\n");
     return TSelectorDraw::Notify();
}

void CMS2Selector::ProcessFill (Long64_t i)
{
//      printf("ProcessFill %lld\n", i);
     cms2.GetEntry(i);
     TSelectorDraw::ProcessFill(i);
}

void CMS2Selector::ProcessFillMultiple (Long64_t i)
{
//      printf("ProcessFillMultiple %lld\n", i);
     cms2.GetEntry(i);
     TSelectorDraw::ProcessFillMultiple(i);
}

void CMS2Selector::ProcessFillObject (Long64_t i)
{
     printf("ProcessFillObject %lld\n", i);
     TSelectorDraw::ProcessFillObject(i);
}

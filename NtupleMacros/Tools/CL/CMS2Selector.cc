// $Id: CMS2Selector.cc,v 1.2 2010/04/18 20:17:03 jmuelmen Exp $

#include <stdio.h>
#include "TChain.h"
#include "CORE/CMS2.h"
#include "CMS2Selector.h"

CMS2Selector::CMS2Selector ()
     : TSelectorDraw(),
       fChain(0)
{

}

void CMS2Selector::Begin (TTree *t)
{
     fChain = dynamic_cast<TChain *>(t);
     iEvent = i_permille_old = 0;
     fEntries = t->GetEntries();
     cms2.Init(t);
     TSelectorDraw::Begin(t);
}

void CMS2Selector::Init (TTree *t)
{
     fChain = dynamic_cast<TChain *>(t);
     TSelectorDraw::Init(t);
     cms2.Init(t);
}

Bool_t CMS2Selector::Notify ()
{
//      printf("Notify\n");
     if (fChain != 0)
	  cms2.Init(fChain->GetTree());
     return TSelectorDraw::Notify();
}

void CMS2Selector::ProcessFill (Long64_t i)
{
//      printf("ProcessFill %lld\n", i);
     int i_permille = (int)floor(1000 * iEvent / float(fEntries));
     if (i_permille != i_permille_old) {
	  i_permille_old = i_permille;
	  // xterm magic from L. Vacavant and A. Cerri
	  if (isatty(1)) {
	       printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
		      "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
	       fflush(stdout);
	  }
     }
     iEvent++;
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

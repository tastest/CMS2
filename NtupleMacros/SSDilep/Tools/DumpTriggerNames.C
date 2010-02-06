// Example that shows how to dump the trigger names
#include <vector>
#include <iostream>
#include <string>
#include "TTree.h"
#include "TBranch.h"
void trigname( TTree* tree )
{
   tree->SetMakeClass(1);   
   std::vector<char> vc;
   if ( TBranch* b = tree->GetBranch("chars_eventMaker_evtHLTtrigNames_CMS2.obj") )
     {
	b->SetAddress(&vc);
	b->GetEntry(0);
	for(unsigned int i=0; i<vc.size(); ++i)
	  std::cout << vc[i];
	std::cout << std::endl;
     }
   else
     cout << "Failed to get the branch" << endl;
}

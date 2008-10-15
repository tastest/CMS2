#include "NMinus1Hist.h"
#include "Sample.h"

using std::vector;
using std::string;

NMinus1Hist::NMinus1Hist (const Sample &s, const std::string &name, 
	     int bins, double min, double max,
	     uint32 cuts_, uint32 cut_mask_)
     : cuts(cuts_),
       h_N(s, name, bins, min, max)
{
     cut_mask.reserve(1);
     mask.reserve(1);
     h_NMinus1.reserve(1);
     cut_mask.push_back(cut_mask_);
     mask.push_back(cuts & ~cut_mask_);
//      printf("cuts: %x, cut_mask: %x, mask: %x\n", cuts, cut_mask[0], mask[0]);
     h_NMinus1[0] = new DileptonHist(s, name + "-N-1", bins, min, max);
}

NMinus1Hist::NMinus1Hist (const Sample &s, const std::string &name, 
	     int bins, double min, double max,
	     uint32 cuts_, vector<uint32> cut_mask_,
	     std::vector<std::string> cut_names)
     : cuts(cuts_),
       cut_mask(cut_mask_),
       h_N(s, name, bins, min, max)
{
     mask.reserve(cut_mask.size());
     h_NMinus1.reserve(cut_mask.size());
     for (vector<uint32>::const_iterator i = cut_mask.begin(); 
	  i != cut_mask.end(); ++i) {
	  mask.push_back(cuts & ~*i);
     }
     for (vector<string>::const_iterator i = cut_names.begin();
	  i != cut_names.end(); ++i) {
	  h_NMinus1.push_back(new DileptonHist(s, name + *i, bins, min, max));
     }
}

void NMinus1Hist::Fill (uint32 cuts_passed, enum DileptonHypType t, double x, 
			double w)
{
     for (unsigned int i = 0; i < mask.size(); ++i) { 
// 	  printf("cuts passed: %x, cuts: %x, cut_mask: %x, mask: %x\t",
// 		 cuts_passed, cuts, cut_mask[i], mask[i]);
	  if ((cuts_passed & mask[i]) == mask[i]) {
// 	       printf("passes N-1\n");
	       h_NMinus1[i]->Fill(t, x, w);
	  } else {
// 	       printf("fails N-1\n");
	  }
     }
     if ((cuts_passed & cuts) == cuts)
	  h_N.Fill(t, x, w);
}

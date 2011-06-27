// Input:  2 histogram
// Output: one histogram which is the efficiency:
// h1 :  TOTAL NUMBER OF EVENTS
// h2 :  NUMBER OF EVENTS THAT PASS

#include "TH1.h"

// Method by pointer
TH1F* eff(TH1F* h1, TH1F* h2, const char* name="eff"){

  // first, verify that all histograms have same binning
  // nx is the number of visible bins
  // nxtot = nx+2 includes underflow and overflow
  Int_t nx = h1->GetNbinsX();
  if (h2->GetNbinsX() != nx) {
    cout << "Histograms must have same number of bins" << endl;
    return 0;
  }

  // get the new histogram
  TH1F* temp = (TH1F*) h1->Clone(name);
  temp->SetTitle(name);
  temp->Reset();
  temp->Sumw2();

  // Do the calculation
  temp->Divide(h2,h1,1.,1.,"B");

  // Done
  return temp;
}


// Method by name
TH1F* eff(const char* name1, const char* name2, const char* name="eff"){

  // Get a list of object and their iterator
  TList* list = gDirectory->GetList() ;
  TIterator* iter = list->MakeIterator();

  // Loop over objects, set the pointers
  TObject* obj;
  TH1F* h1=0;
  TH1F* h2=0;
  TString str1 = Form("%s",name1);
  TString str2 = Form("%s",name2);
  while(obj=iter->Next()) {
    TString objName = obj->GetName();
    if (objName == str1) h1 = (TH1F*) obj;
    if (objName == str2) h2 = (TH1F*) obj;
  }

  // quit if not found
  if (h1 == 0) {
    cout << "Histogram " << name1 << " not found" << endl;
    return 0;
  }
  if (h2 == 0) {
    cout << "Histogram " << name2 << " not found" << endl;
    return 0;
  }

  // Call the method by pointer
  TH1F* temp = eff(h1, h2, name);
  return temp;
}

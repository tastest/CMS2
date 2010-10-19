#ifndef TTBARTHESIS_COMMONFUNCTIONS_C
#define TTBARTHESIS_COMMONFUNCTIONS_C




float GetEntries(const TH1F *h, int lowBin, int highBin) {
  
  if(lowBin < 0 || highBin > h->GetNbinsX()+1) {
    cout << "Bins out of range. Setting the lowBin to the underflow and the highBin to the overflow" << endl;
    lowBin = 0;
    highBin = h->GetNbinsX();
  }
  
  float nentries = 0;
  for(int i = lowBin; i < highBin+1; i++) 
    nentries = nentries + h->GetBinContent(i);
  
  return nentries;
}

float GetTotalError(const TH1F *h, int lowBin, int highBin) {
  
  if(lowBin < 0 || highBin > h->GetNbinsX()+1) {
    cout << "Bins out of range. Setting the lowBin to the underflow and the highBin to the overflow" << endl;
    lowBin = 0;
    highBin = h->GetNbinsX();
  }
  
  float err2 = 0;
  for(int i = lowBin; i < highBin+1; i++) 
    err2 = err2 + pow(h->GetBinError(i),2);
  
  return sqrt(err2);
  
}


std::string formatFloat(double x, const char* formatS) {
  std::string xS = Form(Form("%s", formatS),x);
  double xB = atof(xS.c_str());
  if (x>0 && xB==0){
    xS = Form(" %6.1g",x);
  }
  return xS;
}

#endif

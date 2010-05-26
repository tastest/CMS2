#ifndef ProcDSS_H
#define ProcDSS_H
struct ProcDSChain {
  ProcDSChain(TChain* ch, std::string nm, float sc = 1, bool doW = true, bool chDup = false): 
    events(ch), name(nm), scale1fb(sc),  useWeigtFromBranch(doW),
    checkDuplicates(chDup) {}
  TChain* events;
  std::string name;
  float scale1fb;
  bool useWeigtFromBranch;
  bool checkDuplicates;
};
struct ProcDSS {
  ProcDSS(std::string gnm, float kF = 1, int iPS = 1, bool useCombined = true) : 
    name(gnm), kFactor(kF), prescale(iPS), combine(useCombined) {}
  void add(ProcDSChain& ds) {dsets.push_back(ds);}
  void add(TChain* ch, std::string nm, float sc = 1, bool doW = true, bool chDup = false){
    dsets.push_back(ProcDSChain(ch, nm, sc, doW, chDup));
  }
  unsigned int size() const { return dsets.size();}
  ProcDSChain& operator[](unsigned int iD){
    return dsets[iD];
  }
  std::vector<ProcDSChain> dsets;
  std::string name;
  float kFactor;
  int prescale;
  bool combine;
};

#endif

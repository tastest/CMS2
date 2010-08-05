bool GoodEarlyDataMCTrigger(int dilType){
  bool hlt_ele15_lw_l1r = cms2.passHLTTrigger("HLT_Ele15_LW_L1R");
  bool hltMu9           = cms2.passHLTTrigger("HLT_Mu9");
  // bool hltdiMu3         = cms2.passHLTTrigger("HLT_DoubleMu3");
  //  bool hltdiEle10       = cms2.passHLTTrigger("HLT_DoubleEle10_SWL1R");
  // bool hltdiEle10       = cms2.passHLTTrigger("HLT_DoubleEle5_SW_L1R");

  if (dilType == 0 && ! (hltMu9) ) return false;
  if ((dilType == 1 || dilType == 2) && ! (hltMu9 || hlt_ele15_lw_l1r)) return false;
  if (dilType == 3 && ! hlt_ele15_lw_l1r) return false;

  return true;
}

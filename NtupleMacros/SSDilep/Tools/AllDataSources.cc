
#include "AllDataSources.h"

DataSource fH_DYMM()    { return DataSource("dymm",     H_DYMM, kCyan);    }
DataSource fH_DYEE()    { return DataSource("dyee",     H_DYEE, kMagenta);    }
DataSource fH_DYTT()    { return DataSource("dytt",     H_DYTT, kBlack);    }

DataSource fH_TTBAR()   { return DataSource("ttbar",    H_TTBAR, kYellow);       }
DataSource fH_ZZ()  { return DataSource("zz",   H_ZZ, kGreen);      }
DataSource fH_WZ()  { return DataSource("wz",   H_WZ, kBlue);      }
DataSource fH_WW()      { return DataSource("ww",       H_WW, kRed);          }

DataSource fH_BC20_30()      { return DataSource("QCDBCtoEPt20to30",       H_BC20_30);          }
DataSource fH_BC30_80()      { return DataSource("QCDBCtoEPt30to80",       H_BC30_80);          }
DataSource fH_BC80_170()      { return DataSource("QCDBCtoEPt80to170",       H_BC80_170);          }
DataSource fH_BC30_170()      { return DataSource("QCDBCtoEPt30to170",       H_BC30_170);          }
DataSource fH_BC20_170()      { return DataSource("QCDBCtoEPt20to170",       H_BC20_170);          }

DataSource fH_EM20_30()      { return DataSource("QCDEMenrichedPt20to30",       H_EM20_30);          }
DataSource fH_EM30_80()      { return DataSource("QCDEMenrichedPt30to80",       H_EM30_80);          }
DataSource fH_EM80_170()      { return DataSource("QCDEMenrichedPt80to170",       H_EM80_170);          }
DataSource fH_EM30_170()      { return DataSource("QCDEMenrichedPt30to170",       H_EM30_170);          }
DataSource fH_EM20_170()      { return DataSource("QCDEMenrichedPt20to170",       H_EM20_170);          }

DataSource fH_PHOTONJET() { return DataSource("PhotonJetPt20to170", H_PHOTONJET); }

DataSource fH_WENU()      { return DataSource("we",       H_WENU);          }
DataSource fH_WMUNU()      { return DataSource("wm",       H_WMUNU);          }
DataSource fH_WTAUNU()      { return DataSource("wt",       H_WTAUNU);          }

DataSource fH_QCD30()      { return DataSource("QCDpt30",       H_QCD30);          }
DataSource fH_MU30()       { return DataSource("InclusiveMuPt15",       H_MU30);          };

// 7 TeV V02-00-08 - probably only used by DLE
DataSource fH_WENU_7TeV() { return DataSource("wenu", H_WENU_7TeV); }
DataSource fH_QCD30_7TeV() { return DataSource("qcd_pt30", H_QCD30_7TeV); }
DataSource fH_PHOTONJET_7TeV() { return DataSource("photonjet", H_PHOTONJET_7TeV); }

// old stuff?
DataSource fH_TW()      { return DataSource("tw",       H_TW);          }
DataSource fH_WJETS()   { return DataSource("wjets",    H_WJETS);       }
DataSource fH_MU15_SINGLE() { return DataSource("InclusiveMuPt15", H_MU15_SINGLE, 28); }
DataSource fH_QCD80()      { return DataSource("QCDpt80",       H_QCD80, 17);          }
DataSource fH_WJET_ALP()   { return DataSource("wjetsAlpgen", H_WJET_ALP, 40); }
DataSource fH_ZEEJET_ALP()   { return DataSource("dyeeAlpgen", H_ZEEJET_ALP, 42); }
DataSource fH_ZMMJET_ALP()   { return DataSource("dymmAlpgen", H_ZMMJET_ALP, 44); }
DataSource fH_ZTTJET_ALP()   { return DataSource("dyttAlpgen", H_ZTTJET_ALP, 46); }


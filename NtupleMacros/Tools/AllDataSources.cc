
#include "AllDataSources.h"

DataSource fH_DYMM()    { return DataSource("dymm",     H_DYMM, kAzure-5, "Z/#gamma*#rightarrow #mu^{+}#mu^{-}");    }
DataSource fH_DYEE()    { return DataSource("dyee",     H_DYEE, kAzure-2, "Z/#gamma*#rightarrow e^{+}e^{-}");    }

DataSource fH_DYMM_ODD()    { return DataSource("dymm_odd",     H_DYMM_ODD, kCyan);    }
DataSource fH_DYEE_ODD()    { return DataSource("dyee_odd",     H_DYEE_ODD, kMagenta);    }

DataSource fH_DYMM_EVEN()    { return DataSource("dymm_even",     H_DYMM_EVEN, kCyan);    }
DataSource fH_DYEE_EVEN()    { return DataSource("dyee_even",     H_DYEE_EVEN, kMagenta);    }

DataSource fH_DYTT()    { return DataSource("dytt",     H_DYTT, kAzure+8, "Z/#gamma*#rightarrow#tau^{+}#tau^{-}");    }

DataSource fH_TTBAR()   { return DataSource("ttbar",    H_TTBAR, kRed+1, "t#bar{t}");       }
DataSource fH_ZZ()  { return DataSource("zz",   H_ZZ, kGreen);      }
DataSource fH_WZ()  { return DataSource("wz",   H_WZ, kBlue);      }
DataSource fH_WW()      { return DataSource("ww",       H_WW, kRed);          }
DataSource fH_VV()      { return DataSource("vv",       H_VV, kWhite, "VV");          }

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
DataSource fH_PHOTONJETS() { return DataSource("photonjets", H_PHOTONJETS, kMagenta); }


DataSource fH_WENU()      { return DataSource("we",       H_WENU, kBlue);          }
DataSource fH_WMUNU()      { return DataSource("wm",       H_WMUNU, kBlue + 2);          }
DataSource fH_WTAUNU()      { return DataSource("wt",       H_WTAUNU, kGray);          }

DataSource fH_QCD15()       { return DataSource("QCDpt15",      H_QCD15, kGreen + 2);    };
DataSource fH_QCD30()      { return DataSource("QCDpt30",       H_QCD30, kGreen);          }
DataSource fH_MU30()       { return DataSource("InclusiveMuPt30",       H_MU30);          };
DataSource fH_MU15()       { return DataSource("InclusiveMuPt15",       H_MU15);          };


// old stuff?
DataSource fH_TW()      { return DataSource("tw",       H_TW);          }
DataSource fH_MU15_SINGLE() { return DataSource("InclusiveMuPt15", H_MU15_SINGLE, 28); }
DataSource fH_QCD80()      { return DataSource("QCDpt80",       H_QCD80, kBlue);          }

DataSource fH_WJETS()   { return DataSource("wjets",    H_WJETS,     kGreen-3, "W#rightarrowl#nu");       }
DataSource fH_ZJETS()   { return DataSource("zjets",    H_ZJETS, kCyan);       }
DataSource fH_DATA()      { return DataSource("data",       H_DATA, kRed);          }
DataSource fH_MINBIAS()      { return DataSource("minbias",       H_MINBIAS, kGray);          }

DataSource fH_LM0()   { return DataSource("lm0",    H_LM0, kGray);       }



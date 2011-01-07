
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
DataSource fH_TW()      { return DataSource("tw",       H_TW, kMagenta, "tW");          }
DataSource fH_MU15_SINGLE() { return DataSource("InclusiveMuPt15", H_MU15_SINGLE, 28); }
DataSource fH_QCD80()      { return DataSource("QCDpt80",       H_QCD80, kBlue);          }

DataSource fH_WJETS()   { return DataSource("wjets",    H_WJETS,     kGreen-3, "W#rightarrowl#nu");       }
DataSource fH_ZJETS()   { return DataSource("zjets",    H_ZJETS, kCyan);       }
DataSource fH_DATA()      { return DataSource("data",       H_DATA, kRed);          }
DataSource fH_MINBIAS()      { return DataSource("minbias",       H_MINBIAS, kGray);          }

DataSource fH_LM0()   { return DataSource("lm0",    H_LM0, kGray);       }
DataSource fH_LM1()   { return DataSource("lm1",    H_LM1, kGray);       }
DataSource fH_LM2()   { return DataSource("lm2",    H_LM2, kGray);       }
DataSource fH_LM3()   { return DataSource("lm3",    H_LM3, kGray);       }
DataSource fH_LM4()   { return DataSource("lm4",    H_LM4, kGray);       }
DataSource fH_LM5()   { return DataSource("lm5",    H_LM5, kGray);       }
DataSource fH_LM6()   { return DataSource("lm6",    H_LM6, kGray);       }
DataSource fH_LM7()   { return DataSource("lm7",    H_LM7, kGray);       }
DataSource fH_LM8()   { return DataSource("lm8",    H_LM8, kGray);       }
DataSource fH_LM9()   { return DataSource("lm9",    H_LM9, kGray);       }
DataSource fH_LM10()   { return DataSource("lm10",    H_LM10, kGray);       }

DataSource fH_LQCM300()   { return DataSource("lqcm300",    H_LQCM300, kGray, "LQToCMu_M-300");       }
DataSource fH_LQCM500()   { return DataSource("lqcm500",    H_LQCM500, kGray, "LQToCMu_M-500");       }
DataSource fH_LQUE300()   { return DataSource("lque300",    H_LQUE300, kGray, "LQToUE_M-300");       }
DataSource fH_LQUE500()   { return DataSource("lque500",    H_LQUE500, kGray, "LQToUE_M-500");       }

DataSource fH_D2000LQ1500E()   { return DataSource("d2000lq1500e",    H_D2000LQ1500E, kGray, "D2000LQ1500e");       }
DataSource fH_D1200LQ1100E()   { return DataSource("d1200lq1100e",    H_D1200LQ1100E, kGray, "D1200LQ1100e");       }
DataSource fH_D2000LQ1500M()   { return DataSource("d2000lq1500m",    H_D2000LQ1500M, kGray, "D2000LQ1500m");       }
DataSource fH_D1200LQ1100M()   { return DataSource("d1200lq1100m",    H_D1200LQ1100M, kGray, "D1200LQ1100m");       }

DataSource fH_DQ2500_LDQ2400E()   { return DataSource("dq2500_ldq2400e",    H_DQ2500_LDQ2400E, kGray, "dq2500_ldq2400e");       }
DataSource fH_DQ2500_LDQ2400M()   { return DataSource("dq2500_ldq2400m",    H_DQ2500_LDQ2400M, kGray, "dq2500_ldq2400m");       }
DataSource fH_DQ3500_LDQ1000E()   { return DataSource("dq3500_ldq1000e",    H_DQ3500_LDQ1000E, kGray, "dq3500_ldq1000e");       }
DataSource fH_DQ3500_LDQ1000M()   { return DataSource("dq3500_ldq1000m",    H_DQ3500_LDQ1000M, kGray, "dq3500_ldq1000m");       }



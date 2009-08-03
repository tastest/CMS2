
#include "DataSource.h"

DataSource fH_WW() 	{ return DataSource("ww", 	H_WW); 		}
DataSource fH_TTBAR() 	{ return DataSource("ttbar", 	H_TTBAR); 	}
DataSource fH_DYMM() 	{ return DataSource("dymm", 	H_DYMM);	}
DataSource fH_DYEE() 	{ return DataSource("dyee", 	H_DYEE);	}
DataSource fH_DYTT() 	{ return DataSource("dytt", 	H_DYTT); 	}
DataSource fH_WJETS() 	{ return DataSource("wjets", 	H_WJETS); 	}
DataSource fH_TW()	{ return DataSource("tw", 	H_TW); 		}
DataSource fH_ZZ() 	{ return DataSource("zz", 	H_ZZ); 		}
DataSource fH_WZ() 	{ return DataSource("wz", 	H_WZ); 		}

DataSource fH_EM30_80()      { return DataSource("em30_80",       H_EM30_80);          }
DataSource fH_BC30_80()      { return DataSource("bc30_80",       H_BC30_80);          }
DataSource fH_WENU()      { return DataSource("wenu",       H_WENU);          }

DataSource fH_QCD30()      { return DataSource("QCDpt30",       H_QCD30, 14);          }
DataSource fH_QCD80()      { return DataSource("QCDpt80",       H_QCD80, 17);          }
DataSource fH_WJET_ALP()   { return DataSource("wjetsAlpgen", H_WJET_ALP, 40); }
DataSource fH_ZEEJET_ALP()   { return DataSource("dyeeAlpgen", H_ZEEJET_ALP, 42); }
DataSource fH_ZMMJET_ALP()   { return DataSource("dymmAlpgen", H_ZMMJET_ALP, 44); }
DataSource fH_ZTTJET_ALP()   { return DataSource("dyttAlpgen", H_ZTTJET_ALP, 46); }




#ifndef ALLDATASOURCES_H
#define ALLDATASOURCES_H

#include "DataSource.h"

enum {

        H_WHUNT,

	H_TTBAR   ,
	H_DYMM    ,
	H_DYEE    ,

    H_DYMM_ODD    ,
    H_DYEE_ODD    ,

    H_DYMM_EVEN    ,
    H_DYEE_EVEN    ,


	H_DYTT    ,

	H_ZZ      ,
	H_WZ,
	H_WW      ,

	H_WENU,
	H_WMUNU,
	H_WTAUNU,

	H_EM20_30,
	H_EM30_80,
	H_EM80_170,
	H_EM20_170,
	H_EM30_170,

	H_BC20_30,
	H_BC30_80,
	H_BC80_170,
	H_BC20_170,
	H_BC30_170,

	H_PHOTONJET,
	H_PHOTONJETS,

	H_QCD30,
	H_QCD80,
	H_MU30,
    H_MU15,

	H_TW      ,
	H_WJETS   ,
	H_WJET_ALP,
	H_ZEEJET_ALP,
	H_ZMMJET_ALP,
	H_ZTTJET_ALP,
	H_MU15_SINGLE,

	H_ELEGUNIDEAL,
	H_ELEGUNSTARTUP,

};

/*
const static sources_t sources_all =
(1ll << H_WW)          |
(1ll << H_TTBAR)       |
(1ll << H_DYMM)        |
(1ll << H_DYEE)        |
(1ll << H_DYTT)        |
(1ll << H_WJETS)       |
(1ll << H_TW)          |
(1ll << H_ZZ)          |
(1ll << H_WZ);

const static sources_t sources_dy =
(1ll << H_DYMM)        |
(1ll << H_DYEE);

const static sources_t sources_peaking =
(1ll << H_DYMM)        |
(1ll << H_DYEE)        |
(1ll << H_ZZ)          |
(1ll << H_WZ);

const static sources_t sources_nonpeaking =
(1ll << H_WW)          |
(1ll << H_TTBAR)       |
(1ll << H_DYTT)        |
(1ll << H_WJETS)       |
(1ll << H_TW);
*/

// consistent set of 7 TeV ntuples
// are available with V02-00-12

DataSource fH_WHUNT();

DataSource fH_DYMM();
DataSource fH_DYEE();

DataSource fH_DYMM_ODD();
DataSource fH_DYEE_ODD();

DataSource fH_DYMM_EVEN();
DataSource fH_DYEE_EVEN();

DataSource fH_DYTT();

DataSource fH_TTBAR();
DataSource fH_ZZ();
DataSource fH_WZ();
DataSource fH_WW();

DataSource fH_WENU();
DataSource fH_WMUNU();
DataSource fH_WTAUNU();

DataSource fH_EM20_30();
DataSource fH_EM30_80();
DataSource fH_EM80_170();
DataSource fH_EM30_170();
DataSource fH_EM20_170();

DataSource fH_BC20_30();
DataSource fH_BC30_80();
DataSource fH_BC80_170();
DataSource fH_BC30_170();
DataSource fH_BC20_170();

DataSource fH_PHOTONJET();
DataSource fH_PHOTONJETS();

DataSource fH_QCD30();
DataSource fH_MU30();
DataSource fH_MU15();

// old?
DataSource fH_WJETS();
DataSource fH_TW();
DataSource fH_QCD80();
DataSource fH_WJET_ALP();
DataSource fH_ZEEJET_ALP();
DataSource fH_ZMMJET_ALP();
DataSource fH_ZTTJET_ALP();
DataSource fH_MU15_SINGLE();


DataSource fH_ELEGUNIDEAL();
DataSource fH_ELEGUNSTARTUP();

#endif

